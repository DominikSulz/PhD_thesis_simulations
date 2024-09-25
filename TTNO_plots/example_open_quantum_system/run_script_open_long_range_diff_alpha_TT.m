clear all; clc; close all;

addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\TTNO')
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\TTNO\HSS_and_no_structure_case')
addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\hm-toolbox-master')
addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\tensor_toolbox-master')

%% initializations
% parameters of the model
d = 2^8;                % number of particles
l = log(d)/log(2);      % number of layers
Omega = 0.4;
Delta = -2;
gamma = 1;
alpha = 1;
c_alpha=sum((1:1:d).^(-alpha));
nu = 2/c_alpha;

r_vec = [10 10 10 10 10 10 10 6 4];

% initialisations
sx=[0,1;1,0];     %% Pauli Matrix x
sy=[0,-1i;1i,0];  %% Pauli Matrix y
sz=[1,0;0,-1];    %% Pauli Matrix z
n=[1,0;0,0];      %% Projector onto the excited state Pu=(sz+id)/2;
id=[1,0;0,1];     %% Identity for the single spin
J = [0,0;1,0];

[X,tau] = init_spin_all_dim_diff_rank_TT(d); % initial data - needed for construct exact operator

alpha_save = [0 0.5 1 1.5 2 2.5 3 4 5 6 7 8 9];

hss_tol = [10^-4 10^-12];

for ll=1:length(hss_tol)
    hssoption('threshold',hss_tol(ll));
    for kk=1:length(alpha_save)
        
        alpha = alpha_save(kk);
        
        %% interaction matrix - single particle
        A_single = cell(1,d);
        for ii=1:d
            A_single{ii} = Omega*(-1i*kron(sx,id) + 1i*kron(id,sx.')) ...
                + Delta*(-1i*kron(n,id) + 1i*kron(id,n.')) ...
                + gamma*(kron(J,conj(J)) - 0.5*kron(J'*J,id) - 0.5*kron(id,(J'*J).'));
        end
        V_single = eye(d,d);
        
        %% interaction matrix - long-range interactions
        V_int1 = zeros(d,d);
        A_int1 = cell(1,d);
        V_int2 = zeros(d,d);
        A_int2 = cell(1,d);
        for ii=1:d
            for jj=1:d
                if ii ~= jj
                    V_int1(ii,jj) = -1i*nu * 1/abs(ii-jj)^alpha; % -1i*(nu/2) * 1/abs(ii-jj)^alpha;
                    V_int2(ii,jj) = 1i*nu * 1/abs(ii-jj)^alpha; % 1i*(nu/2) * 1/abs(ii-jj)^alpha;
                end
            end
            A_int1{ii} = kron(n,id);
            A_int2{ii} = kron(id,n.');
        end
        
        %% HSS decompositions
        %% TTNO in TT representation
        cluster = hss();
        cluster = cluster_rec_HSS(cluster,d);
        cluster.topnode = 1;
        
        H_single = hss(V_single,'cluster',cluster);
        H_single = adjust_H(H_single);
        
        H_int1 = hss(V_int1,'cluster',cluster);
        H_int1 = adjust_H(H_int1);
        H_int2 = hss(V_int2,'cluster',cluster);
        H_int2 = adjust_H(H_int2);
        
        %% construction of TTNO with help of HSS
        TTNO_single = TTNO_HSS_construct_single_TT(H_single,A_single,l,l,4*ones(1,d),1:d);
        TTNO_int1 = TTNO_HSS_construct_TT(H_int1,A_int1,l,l,4*ones(1,d),1:d);
        TTNO_int2 = TTNO_HSS_construct_TT(H_int2,A_int2,l,l,4*ones(1,d),1:d);
        
        TTNO = Add_operator_nonround(TTNO_single,TTNO_int1,tau);
        TTNO = Add_operator_nonround(TTNO,TTNO_int2,tau);
        
        TTNO = rounding(TTNO,tau);
        
        %% exact TTNO
        %     B = linearisation_long_range_single(d,J,sx,n,nu,Delta,Omega,gamma,alpha);
        %     TTNO_exact_single = make_operator(X,B,tau,4*ones(1,d));
        %     TTNO_exact_single = rounding(TTNO_exact_single,tau);
        TTNO_exact_single = TTNO_HSS_construct_single_TT(H_single,A_single,l,l,4*ones(1,d),1:d);
        
        TTNO_exact_int1 = TTNO_HSS_construct_TT(A_int1,V_int1,l,l,4*ones(d,1),1:d);
        TTNO_exact_int2 = TTNO_HSS_construct_TT(A_int2,V_int2,l,l,4*ones(d,1),1:d);
        
        TTNO_exact = Add_operator_nonround(TTNO_exact_single,TTNO_exact_int1,tau);
        TTNO_exact = Add_operator_nonround(TTNO_exact,TTNO_exact_int2,tau);
        
        TTNO_exact = rounding(TTNO_exact,tau);
        
        
        %% error check
        tmp = TTNO;
        tmp{end} = -tmp{end};
        E = Add_TTN(TTNO_exact,tmp,tau);
        err(ll,kk) = sqrt(abs(Mat0Mat0(E,E)));
        
        err_scaled(ll,kk) = err(ll,kk)/sqrt(abs(Mat0Mat0(TTNO_exact,TTNO_exact)));
        
        max_rk(ll,kk) = max_rank(TTNO);
        
        ex_rk(ll,kk) = hssrank(H_int1) + hssrank(H_int2) + 4;
        % 4 = 2 + 1 + 1, first 2 comes from H_int1, which contains the
        % identity. First 1 comes from H_int2, where we don't have an identity.
        % Second 1 comes from the diagonal part
        
        kk
    end
end

subplot(1,2,1)
plot(alpha_save,max_rk(1,:),'Linewidth',2)
hold on
plot(alpha_save,ex_rk(1,:),'--','Linewidth',2)
axis([0 alpha_save(end) min(max_rk(1,:))-1 max(max_rk(1,:))+1])
xlabel('{\alpha}','Fontsize',12)
legend('Maximal rank of the TTNO','hssrank \beta_1 + 2 + hssrank \beta_2 + 1 + 1','Fontsize',12)

subplot(1,2,2)
plot(alpha_save,max_rk(2,:),'Linewidth',2)
hold on
plot(alpha_save,ex_rk(2,:),'--','Linewidth',2)
axis([0 alpha_save(end) min(max_rk(2,:))-1 max(max_rk(2,:))+1])
xlabel('{\alpha}','Fontsize',12)
legend('Maximal rank of the TTNO','hssrank \beta_1 + 2 + hssrank \beta_2 + 1 + 1','Fontsize',12)

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