clear all; close all; clc

addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\TTNO')
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\TTNO\HSS_and_no_structure_case')
addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\hm-toolbox-master')
addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\tensor_toolbox-master')


% parameters of the model
Omega = 3;
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

rk = [1 2 3 4 5 6 7 8 9];

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
    
    %% HSS decompositions
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
    
    TTNO = Add_operator_nonround(TTNO_single,TTNO_int1,tau);
    TTNO = Add_operator_nonround(TTNO,TTNO_int2,tau);
    
    TTNO = rounding(TTNO,tau);

    %% exact TTNO
%     B = linearisation_long_range_single(d,J,sx,n,nu,Delta,Omega,gamma,alpha);
%     TTNO_exact_single = make_operator(X,B,tau,4*ones(1,d)); 
%     TTNO_exact_single = rounding(TTNO_exact_single,tau);
    TTNO_exact_single = TTNO_construct_single(A_single,l,l,4*ones(1,d),1:d);
    
    TTNO_exact_int1 = TTNO_no_structure(A_int1,V_int1,l,l,4*ones(d,1),1:d);
    TTNO_exact_int2 = TTNO_no_structure(A_int2,V_int2,l,l,4*ones(d,1),1:d);

    TTNO_exact = Add_operator_nonround(TTNO_exact_single,TTNO_exact_int1,tau);
    TTNO_exact = Add_operator_nonround(TTNO_exact,TTNO_exact_int2,tau);
    
    TTNO_exact = rounding(TTNO_exact,tau);

    
    %% error check
    tmp = TTNO;
    tmp{end} = -tmp{end};
    E = Add_TTN(TTNO_exact,tmp,tau);
%     err(kk) = sqrt(abs(Mat0Mat0(E,E)));
    
    % err_scaled(kk) = err(kk)/sqrt(abs(Mat0Mat0(TTNO_exact,TTNO_exact)));
    
    if kk==9 % to avoid overflow
        E{end} = 10^(-50)*E{end};
        err(kk) = sqrt(abs(Mat0Mat0(E,E)));
        TTNO_exact{end} = 10^(-50)*TTNO_exact{end};
        err_scaled(kk) = err(kk)/sqrt(abs(Mat0Mat0(TTNO_exact,TTNO_exact)));
    else
        err(kk) = sqrt(abs(Mat0Mat0(E,E)));
        err_scaled(kk) = err(kk)/sqrt(abs(Mat0Mat0(TTNO_exact,TTNO_exact)));
    end
    
    max_rk(kk) = max_rank(TTNO);
    
    ex_rk(kk) = hssrank(H_int1) + hssrank(H_int2) + 4;  
    % 4 = 2 + 1 + 1, first 2 comes from H_int1, which contains the
    % identity. First 1 comes from H_int2, where we don't have an identity.
    % Second 1 comes from the diagonal part
    
    kk
    
end

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
legend('Maximal tree rank','Expected maximal tree rank','Fontsize',20)



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


% %%% check individual error
% tmp = TTNO_single;
% tmp{end} = -tmp{end};
% E = Add_TTN(tmp,TTNO_exact_single,tau);
% sqrt(Mat0Mat0(E,E))/sqrt(Mat0Mat0(TTNO_exact_single,TTNO_exact_single))
% 
% tmp = TTNO_int1;
% tmp{end} = -tmp{end};
% E = Add_TTN(tmp,TTNO_exact_int1,tau);
% err_int1 = sqrt(Mat0Mat0(E,E))/sqrt(Mat0Mat0(TTNO_exact_int1,TTNO_exact_int1))
% 
% tmp = TTNO_int2;
% tmp{end} = -tmp{end};
% E = Add_TTN(tmp,TTNO_exact_int2,tau);
% err_int2 = sqrt(Mat0Mat0(E,E))/sqrt(Mat0Mat0(TTNO_exact_int2,TTNO_exact_int2))
% 
% T1 = Add_TTN(TTNO_int1,TTNO_int2,tau);
% T2 = Add_TTN(TTNO_exact_int1,TTNO_exact_int2,tau);
% tmp = T1;
% tmp{end} = -tmp{end};
% E = Add_TTN(tmp,T2,tau);
% err_int12 = sqrt(Mat0Mat0(E,E))/sqrt(Mat0Mat0(T2,T2))
% 
% T1 = Add_TTN(TTNO_int1,TTNO_single,tau);
% T2 = Add_TTN(TTNO_exact_int1,TTNO_exact_single,tau);
% tmp = T1;
% tmp{end} = -tmp{end};
% E = Add_TTN(tmp,T2,tau);
% err_int1s = sqrt(Mat0Mat0(E,E))/sqrt(Mat0Mat0(T2,T2))
% 
% T1 = Add_TTN(TTNO_int2,TTNO_single,tau);
% T2 = Add_TTN(TTNO_exact_int2,TTNO_exact_single,tau);
% tmp = T1;
% tmp{end} = -tmp{end};
% E = Add_TTN(tmp,T2,tau);
% err_int2s = sqrt(Mat0Mat0(E,E))/sqrt(Mat0Mat0(T2,T2))
% 
% T1 = Add_TTN(TTNO_int1,TTNO_single,tau);
% T1 = Add_TTN(T1,TTNO_int2,tau);
% T2 = Add_TTN(TTNO_exact_int1,TTNO_exact_single,tau);
% T2 = Add_TTN(T2,TTNO_exact_int2,tau);
% tmp = T1;
% tmp{end} = -tmp{end};
% E = Add_TTN(tmp,T2,tau);
% err_int12s = sqrt(Mat0Mat0(E,E))/sqrt(Mat0Mat0(T2,T2))

%%%
% T1 = Add_operator_nonround(TTNO_int1,TTNO_single,tau);
% T1 = Add_operator_nonround(T1,TTNO_int2,tau);
% T2 = Add_operator_nonround(TTNO_exact_int1,TTNO_exact_single,tau);
% T2 = Add_operator_nonround(T2,TTNO_exact_int2,tau);
% tmp = T1;
% tmp{end} = -tmp{end};
% E = Add_TTN(tmp,T2,tau);
% err_int12s_nr = sqrt(Mat0Mat0(E,E))/sqrt(Mat0Mat0(T2,T2))
% 
% T1 = Add_operator_nonround(TTNO_int1,TTNO_single,tau);
% T1 = rounding(T1,tau);
% T1 = Add_operator_nonround(T1,TTNO_int2,tau);
% T1 = rounding(T1,tau);
% T2 = Add_operator_nonround(TTNO_exact_int1,TTNO_exact_single,tau);
% T2 = rounding(T2,tau);
% T2 = Add_operator_nonround(T2,TTNO_exact_int2,tau);
% T2 = rounding(T2,tau);
% tmp = T1;
% tmp{end} = -tmp{end};
% E = Add_TTN(tmp,T2,tau);
% err_int12s_nr2 = sqrt(Mat0Mat0(E,E))/sqrt(Mat0Mat0(T2,T2))

% %%%