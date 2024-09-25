clear all; close all; clc

addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\TTNO')
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\TTNO\HSS_and_no_structure_case')
addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\hm-toolbox-master')
addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\tensor_toolbox-master')

% addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\rank_adaptive_integrator_for_TTN')

% parameters of the model
Omega = 3;
Delta = -2;
gamma = 1;
alpha = 1;

r_vec = [4 4 4 4 4 4 4 4 4]; 

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
    nu = 2/c_alpha; % 2/c_alpha; % 0.5/c_alpha; paper with Federico
    
    [X,tau] = init_diss_all_dim_diff_rank(r_vec,2,d); % initial data - binary balanced tree
%     X = truncate(X,10^-12,4,4); % make it rank 4
    [X_TT,tau_TT] = init_spin_all_dim_diff_rank_TT(d); % initial data - MPS
    
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

    
    %% HSS
    H_single = hss(V_single,'cluster',1:d);
    H_single = adjust_H(H_single);
    
    H_int1 = hss(V_int1,'cluster',1:d);
    H_int1 = adjust_H(H_int1);
    H_int2 = hss(V_int2,'cluster',1:d);
    H_int2 = adjust_H(H_int2);
    
    %% construction of TTNO with help of HSS
    TTNO_single = TTNO_HSS_construct_single(H_single,A_single,l,l,4*ones(1,d),1:d);
%     TTNO_single = rounding(TTNO_single,tau);
    TTNO_int1 = TTNO_HSS_construct(H_int1,A_int1,l,l,4*ones(1,d),1:d);
%     TTNO_int1 = rounding(TTNO_int1,tau);
    TTNO_int2 = TTNO_HSS_construct(H_int2,A_int2,l,l,4*ones(1,d),1:d);
%     TTNO_int2 = rounding(TTNO_int2,tau);
    
    TTNO = Add_operator_nonround(TTNO_single,TTNO_int1,tau);
%     TTNO = rounding(TTNO,tau);
    TTNO = Add_operator_nonround(TTNO,TTNO_int2,tau);
    
%     TTNO = Add_TTN(TTNO_single,TTNO_int1,tau);
%     TTNO = Add_TTN(TTNO,TTNO_int2,tau);
    
    
%     tmp = rounding(TTNO,tau);
%     tmp{end} = - tmp{end};
%     E = Add_operator_nonround(TTNO,tmp,tau);
%     err = Mat0Mat0(E,E)
%     err = Mat0Mat0(E,E)/Mat0Mat0(TTNO,TTNO)
%
% %     %%%
%     test = TTNO;
%     test{end} = -test{end};
% %     TTNO{end} = (1/Mat0Mat0(TTNO,TTNO))*TTNO{end};
    TTNO = rounding(TTNO,tau);
% %     TTNO = rounding2(TTNO,tau);
% %       TTNO = rounding_root2leaf(TTNO,tau);
%     TTNO = truncate(TTNO,10^-14,200,2);
%     %%%
%     E = Add_operator_nonround(TTNO,test,tau);
%     Mat0Mat0(E,E)
%     Mat0Mat0(E,E)/Mat0Mat0(TTNO,TTNO)

%% TTNO in TT representation
    cluster = hss();
    cluster = cluster_rec_HSS(cluster,d);
    cluster.topnode = 1;
    
    H_single_TT = hss(V_single,'cluster',cluster);
    H_int1_TT = hss(V_int1,'cluster',cluster);
    H_int1_TT = adjust_H(H_int1_TT);
    H_int2_TT = hss(V_int2,'cluster',cluster);
    H_int2_TT = adjust_H(H_int2_TT);
    
    TTNO_TT_single = TTNO_HSS_construct_single_TT(H_single_TT,A_single,d-1,d-1,4*ones(1,d),1:d);
%     TTNO_TT_single = rounding(TTNO_TT_single,tau_TT);
    TTNO_TT_int1 = TTNO_HSS_construct_TT(H_int1_TT,A_int1,d-1,d-1,4*ones(1,d),1:d);
%     TTNO_TT_int1 = rounding(TTNO_TT_int1,tau_TT);
    TTNO_TT_int2 = TTNO_HSS_construct_TT(H_int2_TT,A_int2,d-1,d-1,4*ones(1,d),1:d);
%     TTNO_TT_int2 = rounding(TTNO_TT_int2,tau_TT);
    
    TTNO_TT = Add_operator_nonround(TTNO_TT_single,TTNO_TT_int1,tau_TT);
%     TTNO_TT = rounding2(TTNO_TT,tau_TT);
    TTNO_TT = Add_operator_nonround(TTNO_TT,TTNO_TT_int2,tau_TT);
    
%     TTNO_TT = Add_TTN(TTNO_TT_single,TTNO_TT_int1,tau_TT);
% %     TTNO_TT = rounding2(TTNO_TT,tau_TT);
%     TTNO_TT = Add_TTN(TTNO_TT,TTNO_TT_int2,tau_TT);
    
%     tmp_TT = TTNO_TT;
    
%     tmp_TT = rounding(TTNO_TT,tau_TT);
%     tmp_TT{end} = - tmp_TT{end};
%     E_TT = Add_operator_nonround(TTNO_TT,tmp_TT,tau_TT);
%     err_TT = Mat0Mat0(E_TT,E_TT)
%     err_TT = Mat0Mat0(E_TT,E_TT)/Mat0Mat0(TTNO_TT,TTNO_TT)
   
%
% % %     %%%
%     test_TT = TTNO_TT;
%     test_TT{end} = -test_TT{end};
%     TTNO_TT{end} = (1/Mat0Mat0(TTNO_TT,TTNO_TT))*TTNO_TT{end};
    TTNO_TT = rounding(TTNO_TT,tau_TT);
% %     TTNO_TT = rounding2(TTNO_TT,tau_TT);
% %     TTNO_TT = rounding_root2leaf(TTNO_TT,tau_TT);
%     TTNO_TT = truncate(TTNO_TT,10^-14,200,2);
%     %%%
%     E = Add_operator_nonround(TTNO_TT,test_TT,tau_TT);
%     Mat0Mat0(E,E)
%     Mat0Mat0(E,E)/Mat0Mat0(TTNO_TT,TTNO_TT)
    
    %% max ranks
    max_rk(kk) = max_rank(TTNO);
    max_rk_TT(kk) = max_rank(TTNO_TT);
    
    ex_rk(kk) = hssrank(H_int1) + hssrank(H_int2) + 4;  
    % 4 = 2 + 1 + 1, first 2 comes from H_int1, which contains the
    % identity. First 1 comes from H_int2, where we don't have an identity.
    % Second 1 comes from the diagonal part
    ex_rk_TT(kk) = hssrank(H_int1_TT) + hssrank(H_int2_TT) + 4;
    
%     % summend ranks
%     sum_rank(kk) = sum_ranks(TTNO);
%     sum_rank_TT(kk) = sum_ranks(TTNO_TT);
    
    % memory
    S = whos('TTNO');
    mem(kk) = S.bytes;
    S = whos('TTNO_TT');
    mem_TT(kk) = S.bytes;
    
    % application
    tmp1 = apply_operator_nonglobal(X,TTNO,d);
    S = whos('tmp1');
    mem_app(kk) = S.bytes;
    tmp2 = apply_operator_nonglobal(X_TT,TTNO_TT,d);
    S = whos('tmp2');
    mem_app_TT(kk) = S.bytes;
    
%     % free parameters
%     free_p(kk) = count_free_parameter(TTNO);
%     free_p_TT(kk) = count_free_parameter(TTNO_TT);

    % product ranks
    [sum_rkp,N] = prod_ranks(TTNO);
    prod_rk(kk) = sum_rkp/N;
    [sum_rkp_TT,N_TT] = prod_ranks(TTNO_TT);
    prod_rk_TT(kk) = sum_rkp_TT/N_TT;
    
    % summend ranks
    sum_rk(kk) = sum_ranks(TTNO);
    sum_rk_TT(kk) = sum_ranks(TTNO_TT);
    
    kk
    
end

% plot
figure(1)
subplot(1,2,1)
plot(2.^rk,max_rk,'x-','Linewidth',2)
hold on
plot(2.^rk,max_rk_TT,'*--','Linewidth',2)
% plot(2.^rk,ex_rk,'-.','Linewidth',2)
% plot(2.^rk,ex_rk_TT,':','Linewidth',2)
xlabel('Number of particles','Fontsize',20)
legend('Balanced binary tree','Tensor train',...
'Fontsize',20)
title('Maximal tree rank','Fontsize',20)

subplot(1,2,2)
plot(2.^rk,mem,'x-','Linewidth',2)
hold on
plot(2.^rk,mem_TT,'*--','Linewidth',2)
xlabel('Number of particles','Fontsize',20)
legend('Balanced binary tree','Tensor train'...
      ,'Fontsize',20)
  title('Memory in bytes','Fontsize',20)

% figure(2)
% plot(2.^rk,free_p,'Linewidth',2)
% hold on
% plot(2.^rk,free_p_TT,'--','Linewidth',2)
% xlabel('Number of particles','Fontsize',12)
% legend('Free paramters TTNO binary tree','Free paramters TTNO TT'...
%       ,'Fontsize',12)

figure(2)
plot(2.^rk,mem_app,'Linewidth',2)
hold on
plot(2.^rk,mem_app_TT,'--','Linewidth',2)
xlabel('Number of particles','Fontsize',12)
legend('Memory application TTNO binary tree','Memory application TTNO TT'...
      ,'Fontsize',12)
  
figure(3)
plot(2.^rk,prod_rk,'Linewidth',2)
hold on
plot(2.^rk,prod_rk_TT,'--','Linewidth',2)
xlabel('Number of particles','Fontsize',12)
legend('Averaged product ranks TTNO binary tree','Averaged product ranks TTNO TT'...
      ,'Fontsize',12)
  
figure(4)
plot(2.^rk,sum_rk,'Linewidth',2)
hold on
plot(2.^rk,sum_rk_TT,'--','Linewidth',2)
xlabel('Number of particles','Fontsize',12)
legend('Summend ranks TTNO binary tree','Summed ranks TTNO TT'...
      ,'Fontsize',12)

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

% %% error check of MPS
% % build up the full operator
% B = cell(1,d);
% for ii=1:d
%     B{ii,ii} = Omega*(-1i*kron(sx,id) + 1i*kron(id,sx.')) ...
%         + Delta*(-1i*kron(n,id) + 1i*kron(id,n.')) ...
%         + gamma*(kron(J,conj(J)) - 0.5*kron(J'*J,id) - 0.5*kron(id,(J'*J).'));
% end
% count = d + 1;
% for ii=1:d
%     for jj=1:d
%         if ii ~= jj
%             B{count,ii} = 0.5* -1i*nu * 1/abs(ii-jj)^alpha *kron(n,id);
%             B{count,jj} = kron(n,id);
%             count = count + 1;
%         end
%     end
% end
% for ii=1:d
%     for jj=1:d
%         if ii ~= jj
%             B{count,ii} = 0.5* 1i*nu * 1/abs(ii-jj)^alpha * kron(id,n.');
%             B{count,jj} = kron(id,n.');
%             count = count + 1;
%         end
%     end
% end
% H_ex = make_operator(X_TT,B,tau_TT,4*ones(d,1));
% H_ex = rounding(H_ex,tau_TT);
% 
% H_ex{end} = -H_ex{end};
% D = Add_TTN(H_ex,TTNO_TT,tau_TT);
% err = sqrt(abs(Mat0Mat0(D,D)))

% old code for TT construction
% with hss code
% z = sparse([]);
% %     % old
% for jj=1:d-1
%     tmp = jj*ones(1,2^(d-jj-1));
%     z = [z tmp];
% end
% z(end+1) = d;
%
% %     z = sparse([]);
% %     for jj=1:d-1
% %         if (2^(d-jj-1) - 2) > 0
% %             z = [z jj sparse(1,2^(d-jj-1) - 2) jj];
% %         elseif (2^(d-jj-1) - 2) == 0
% %             z = [z jj jj];
% %         else
% %             z = [z jj];
% %         end
% %     end
% %     z = [z d];
%
% %     z = 1:d;



% old code
% % test with own hss code
% %     z = [];
% % %     % old
% %     for jj=1:d
% %         tmp = [jj jj];
% %         z = [z tmp];
% %     end
% z = 1:d;
%
% H_single_TT2 = hss_TT(V_single,'cluster',z);
% H_int1_TT2 = hss_TT(V_int1,'cluster',z);
% H_int2_TT2 = hss_TT(V_int2,'cluster',z);
%
% TTNO_TT_single2 = TTNO_HSS_construct_single_TT(H_single_TT2,A_single,l,l,4*ones(1,d),1:d);
% TTNO_TT_int12 = TTNO_HSS_construct_TT(H_int1_TT2,A_int1,l,l,4*ones(1,d),1:d);
% TTNO_TT_int22 = TTNO_HSS_construct_TT(H_int2_TT2,A_int2,l,l,4*ones(1,d),1:d);
%
% TTNO_TT2 = Add_TTN(TTNO_TT_single2,TTNO_TT_int12,tau_TT);
% %     TTNO_TT = rounding(TTNO_TT,tau_TT);
% TTNO_TT2 = Add_TTN(TTNO_TT2,TTNO_TT_int22,tau_TT);
% TTNO_TT2 = rounding(TTNO_TT2,tau_TT);
% %     TTNO_TT = truncate(TTNO_TT,10^-10,100,2);