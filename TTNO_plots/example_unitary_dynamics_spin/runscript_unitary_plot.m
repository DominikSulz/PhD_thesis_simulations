clear all; clc; close all;

addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\TTNO')
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\TTNO\HSS_and_no_structure_case')
addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\hm-toolbox-master')
addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\tensor_toolbox-master')

%% initializations
n = 2;             % physical dimension

sx=[0,1;1,0];      % Pauli Matrix \sigma_x
n_Pu=[1,0;0,0];    % Projector onto the excited state Pu=(sz+id)/2;
nu = 2;
Delta = -2;
Omega = 3;
alpha = 1;

rk = [3 4 5 6 7 8 9];

for kk=1:length(rk)
    d = 2^rk(kk);           % number of particles
    l = log(d)/log(2);      % number of layers
    
    [X,tau] = init_spin_all_dim_same_rank(4,1,d); % initial data - needed for construct exact operator
    
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
    hssoption('threshold',10^-12);
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
    TTNO_exact_single = rounding(TTNO_exact_single,tau);
    
    TTNO_exact_int = TTNO_no_structure(A_int,V_int,l,l,n*ones(d,1),1:d);
    TTNO_exact_int = rounding(TTNO_exact_int,tau);
    
    TTNO_exact = Add_TTN(TTNO_exact_single,TTNO_exact_int,tau);
    
    
    %% error check
    tmp = TTNO;
    tmp{end} = -tmp{end};
    E = Add_TTN(TTNO_exact,tmp,tau);
    err(kk) = sqrt(abs(Mat0Mat0(E,E)));
    
    err_scaled(kk) = sqrt(abs(Mat0Mat0(E,E)))/sqrt(abs(Mat0Mat0(TTNO_exact,TTNO_exact)));
    
    max_rk(kk) = max_rank(TTNO);
    
    ex_rk(kk) = hssrank(H_int) + 2 + 1; % +1 for the non interacting part
    
    
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
hssoption('threshold',10^-12);

% plot
subplot(1,2,1)
semilogy(2.^rk,err_scaled,'x-','Linewidth',2)
xlabel('Number of particles','Fontsize',20)
legend('Relative error','Fontsize',20)

subplot(1,2,2)
plot(2.^rk,max_rk,'x-','Linewidth',2)
hold on
plot(2.^rk,ex_rk,'*--','Linewidth',2)
xlabel('Number of particles','Fontsize',20)
legend('Maximal tree rank','HSS rank + 3','Fontsize',20)
% legend('Maximal rank of the TTNO')