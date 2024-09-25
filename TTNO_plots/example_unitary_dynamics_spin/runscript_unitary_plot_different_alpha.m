clear all; clc; close all;

addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\TTNO')
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\TTNO\HSS_and_no_structure_case')
addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\hm-toolbox-master')
addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\tensor_toolbox-master')

%% initializations
n = 2;             % physical dimension
d = 2^8;           % number of particles
l = log(d)/log(2); % number of layers

sx=[0,1;1,0];      % Pauli Matrix \sigma_x
n_Pu=[1,0;0,0];    % Projector onto the excited state Pu=(sz+id)/2;
nu = 2;
Delta = -2;
Omega = 3;

[X,tau] = init_spin_all_dim_same_rank(4,1,d); % initial data - needed for construct exact operator

alpha_save = [0 0.5 1 1.5 2 2.5 3 4 5 6 7 8 9];
hss_tol = [10^-6 10^-12];

for ll=1:length(hss_tol)
    hssoption('threshold',hss_tol(ll));
    for kk=1:length(alpha_save)
        
        alpha = alpha_save(kk);
        
        %% interaction matrix - single particle
        A_single = cell(1,d);
        for ii=1:d
            A_single{ii} = Omega*sx + Delta*n_Pu;
        end
        V_single = eye(d,d);
        
        %% interaction matrix - long-range interactions
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
        
        %% HSS
        H_single = hss(V_single,'cluster',1:d);
        H_single = adjust_H(H_single);
        
        H_int = hss(V_int,'cluster',1:d);
        H_int = adjust_H(H_int);
        
        %% construction of TTNO with help of HSS
        TTNO_single = TTNO_HSS_construct_single(H_single,A_single,l,l,n*ones(1,d),1:d);
        TTNO_int = TTNO_HSS_construct(H_int,A_int,l,l,n*ones(1,d),1:d);
        
        TTNO = Add_TTN(TTNO_single,TTNO_int,tau);
        TTNO = rounding(TTNO,tau);
        
        %% exact TTNO
        B = linearisation_long_range_unitary(d,sx,n_Pu,nu,Delta,Omega,alpha);
        TTNO_exact_single = make_operator(X,B,tau,n*ones(1,d));
        
        TTNO_exact_int = TTNO_no_structure(A_int,V_int,l,l,n*ones(d,1),1:d);
        
        TTNO_exact = Add_TTN(TTNO_exact_single,TTNO_exact_int,tau);
        
        
        %% error check
        tmp = TTNO;
        tmp{end} = -tmp{end};
        E = Add_TTN(TTNO_exact,tmp,tau);
        err(ll,kk) = sqrt(abs(Mat0Mat0(E,E)));
        
        err_scaled(ll,kk) = sqrt(abs(Mat0Mat0(E,E)))/sqrt(abs(Mat0Mat0(TTNO_exact,TTNO_exact)));
        
        max_rk(ll,kk) = max_rank(TTNO);
        
        ex_rk(ll,kk) = hssrank(H_int) + 2 + 1; % +1 for the non interacting part
        
        
        %     % extra check
        %     B = linearisation_long_range_unitary_full(d,sx,n_Pu,nu,Delta,Omega,alpha);
        %     TTNO_exact_test = make_operator(X,B,tau,n*ones(1,d));
        %     TTNO_exact_test = rounding(TTNO_exact_test,tau);
        %     tmp = TTNO_exact;
        %     tmp{end} = -tmp{end};
        %     E = Add_TTN(TTNO_exact_test,tmp,tau);
        %     sqrt(abs(Mat0Mat0(E,E)))
        
        kk
    end
end

% plot
subplot(1,2,1)
plot(alpha_save,max_rk(1,:),'x-','Linewidth',2)
hold on
plot(alpha_save,ex_rk(1,:),'*--','Linewidth',2)
axis([0 alpha_save(end) min(max_rk(1,:))-1 max(max_rk(1,:))+1])
xlabel('{\alpha}','Fontsize',20)
legend('Maximal tree rank','hssrank + 3','Fontsize',20)

subplot(1,2,2)
plot(alpha_save,max_rk(2,:),'x-','Linewidth',2)
hold on
plot(alpha_save,ex_rk(2,:),'*--','Linewidth',2)
axis([0 alpha_save(end) min(max_rk(2,:))-1 max(max_rk(2,:))+1])
xlabel('{\alpha}','Fontsize',20)
legend('Maximal tree rank','hssrank + 3','Fontsize',20)