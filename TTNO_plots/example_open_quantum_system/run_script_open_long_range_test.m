clear all; close all; clc

addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\TTNO')
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\TTNO\HSS_and_no_structure_case')
addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\hm-toolbox-master')
addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\tensor_toolbox-master')

addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\rank_adaptive_integrator_for_TTN')

% the code works now and gives good results. 2 things are to be discussed:
% 1) Depending on the parameters, the Frobeniusnorm of the TTNO get super
% big (10^80) such that svd etc. have trouble and potentially crash. 
% 2) Probably linked to 1). The error in scaled Frobeniusnorm is not of maschine
% precision anymore when taking many particles. Probably because the norm
% in total is to big.

% parameters of the model
Omega = 0.4;
Delta = -2;
gamma = 1; 
alpha = 1;

r_vec = [10 10 10 10 10 10 10 6 4]; 

% initialisations
sx=[0,1;1,0];     %% Pauli Matrix x
sy=[0,-1i;1i,0];  %% Pauli Matrix y
sz=[1,0;0,-1];    %% Pauli Matrix z
n=[1,0;0,0];      %% Projector onto the excited state Pu=(sz+id)/2;
id=[1,0;0,1];     %% Identity for the single spin
J = [0,0;1,0];

rk = [1 2 3 4 5 6];

for kk=1:length(rk)
    d = 2^rk(kk);           % number of particles
    l = log(d)/log(2);      % number of layers
    
    c_alpha=sum((1:1:d).^(-alpha));
    nu = 2/c_alpha;
    
    [X,tau] = init_diss_all_dim_diff_rank(r_vec,2,d); % initial data - needed for construct exact operator
    
    % test
    % [X,tau] = init_unconventional(r,2,d);
    
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

    % *0.5 is away, as for the TTNO we only use the upper triangluar part of V,
    % while in the direct construction we use the full matrix V
    % --> also change that in the draft!
    
    %% HSS
    H_single = hss(V_single,'cluster',1:d);
    H_single = adjust_H(H_single);
    
    H_int1 = hss(V_int1,'cluster',1:d);
    H_int1 = adjust_H(H_int1);
    H_int2 = hss(V_int2,'cluster',1:d);
    H_int2 = adjust_H(H_int2);
    
    %% construction of TTNO with help of HSS
    TTNO_single = TTNO_HSS_construct_single(H_single,A_single,l,l,4*ones(1,d),1:d);
    TTNO_int1 = TTNO_HSS_construct(H_int1,A_int1,l,l,4*ones(1,d),1:d);
    TTNO_int2 = TTNO_HSS_construct(H_int2,A_int2,l,l,4*ones(1,d),1:d);
    
    TTNO = Add_TTN(TTNO_single,TTNO_int1,tau);
    TTNO = rounding(TTNO,tau);
    TTNO = Add_TTN(TTNO,TTNO_int2,tau);
    TTNO = rounding(TTNO,tau);
   
%     TTNO = truncate(TTNO,10^-14,100,2); % eliminates one of the identities, 
                                          % as we have one from H_int1 and one from H_int2
    
    %% exact TTNO
    B = linearisation_long_range_single(d,J,sx,n,nu,Delta,Omega,gamma,alpha);
    TTNO_exact_single = make_operator(X,B,tau,4*ones(1,d));
    TTNO_exact_single = rounding(TTNO_exact_single,tau);
    
    TTNO_exact_int1 = TTNO_no_structure(A_int1,V_int1,l,l,4*ones(d,1),1:d);
%     TTNO_exact_int1 = rounding(TTNO_exact_int1,tau);
    TTNO_exact_int2 = TTNO_no_structure(A_int2,V_int2,l,l,4*ones(d,1),1:d);
%     TTNO_exact_int2 = rounding(TTNO_exact_int2,tau);
    
    TTNO_exact = Add_TTN(TTNO_exact_single,TTNO_exact_int1,tau);
    TTNO_exact = rounding(TTNO_exact,tau);
    TTNO_exact = Add_TTN(TTNO_exact,TTNO_exact_int2,tau);
    TTNO_exact = rounding(TTNO_exact,tau);
    
%      --> Addition von beiden long range teilen ist Problem! 
    %% error check
    tmp = TTNO;
    tmp{end} = -tmp{end};
    E = Add_TTN(TTNO_exact,tmp,tau);
%     E = rounding(E,tau);
    err(kk) = sqrt(abs(Mat0Mat0(E,E)));
    
    err_scaled(kk) = err(kk)/sqrt(abs(Mat0Mat0(TTNO_exact,TTNO_exact)));
    
    max_rk(kk) = max_rank(TTNO);
    
    ex_rk(kk) = hssrank(H_int1) + 2 + hssrank(H_int2) + 2 + 1;
    
    
    kk
    
% %     % extra check
%     B = linearisation_long_range(d,J,sx,n,nu,Delta,Omega,gamma,alpha);
%     TTNO_exact_test = make_operator(X,B,tau,4*ones(1,d));
%     TTNO_exact_test = rounding(TTNO_exact_test,tau);
%     tmp = TTNO_exact;
%     tmp{end} = -tmp{end};
%     E = Add_TTN(TTNO_exact_test,tmp,tau);
%     sqrt(abs(Mat0Mat0(E,E)))/sqrt(abs(Mat0Mat0(TTNO_exact_test,TTNO_exact_test)))
    
%     B = linearisation_long_range_test(d,J,sx,n,nu,Delta,Omega,gamma,alpha);
%     TTNO_exact_int1_test = make_operator(X,B,tau,4*ones(1,d));
%     TTNO_exact_int1_test{end} = 0.5*TTNO_exact_int1_test{end}; 
    
end

% plot
subplot(1,2,1)
semilogy(2.^rk,err_scaled,'Linewidth',2)
xlabel('Number of particles','Fontsize',12)
legend('Scaled error in Frobenius norm','Fontsize',12)

subplot(1,2,2)
plot(2.^rk,max_rk,'Linewidth',2)
hold on
plot(2.^rk,ex_rk,'--','Linewidth',2)
xlabel('Number of particles','Fontsize',12)
legend('Maximal rank of the TTNO','hssrank + 4 + 1','Fontsize',12)
