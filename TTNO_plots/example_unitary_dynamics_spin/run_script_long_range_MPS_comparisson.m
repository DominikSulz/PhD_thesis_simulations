clear all; close all; clc

addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\TTNO')
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\TTNO\HSS_and_no_structure_case')
addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\hm-toolbox-master')
addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\tensor_toolbox-master')

% addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\rank_adaptive_integrator_for_TTN')

% parameters of the model
n = 2;             % physical dimension

Delta = -2;
Omega = 3;
alpha = 1;
nu = 2;

r_vec = [4 4 4 4 4 4 4 4 4];

% initialisations
sx=[0,1;1,0];     %% Pauli Matrix x
sy=[0,-1i;1i,0];  %% Pauli Matrix y
sz=[1,0;0,-1];    %% Pauli Matrix z
n_Pu=[1,0;0,0];    % Projector onto the excited state Pu=(sz+id)/2;
id=[1,0;0,1];     %% Identity for the single spin
J = [0,0;1,0];

rk = [1 2 3 4 5 6 7 8 9];

for kk=1:length(rk)
    d = 2^rk(kk);           % number of particles
    l = log(d)/log(2);      % number of layers
    
    [X,tau] = init_spin_all_dim_same_rank(r_vec,2,d); % initial data - binary balanced tree
    [X_TT,tau_TT] = init_spin_all_dim_diff_rank_TT(d); % initial data - MPS
    
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
    
    %% TTNO in TT representation
    cluster = hss();
    cluster = cluster_rec_HSS(cluster,d);
    cluster.topnode = 1;
    
    H_single_TT = hss(V_single,'cluster',cluster);
    H_int_TT = hss(V_int,'cluster',cluster);
    
    TTNO_TT_single = TTNO_HSS_construct_single_TT(H_single_TT,A_single,d-1,d-1,n*ones(1,d),1:d);
    TTNO_TT_int = TTNO_HSS_construct_TT(H_int_TT,A_int,d-1,d-1,n*ones(1,d),1:d);
    
    TTNO_TT = Add_TTN(TTNO_TT_single,TTNO_TT_int,tau_TT);
    TTNO_TT = rounding(TTNO_TT,tau_TT);
    %     TTNO_TT = rounding(TTNO_TT,tau_TT);
    %     TTNO_TT = truncate(TTNO_TT,10^-10,100,2);
    
    %% max ranks
    max_rk(kk) = max_rank(TTNO);
    max_rk_TT(kk) = max_rank(TTNO_TT);
    
%     % predicted ranks
%     max_rk_ex(kk) = hssrank(H_int) + 3;
%     max_rk_TT_ex(kk) = hssrank(H_int_TT) + 3;
%     
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
    
    kk
    
end

% plot
figure(1)
subplot(1,2,1)
plot(2.^rk,max_rk,'x-','Linewidth',2)
hold on
plot(2.^rk,max_rk_TT,'*--','Linewidth',2)
% plot(2.^rk,max_rk_ex,'--','Linewidth',2)
% plot(2.^rk,max_rk_TT_ex,'--','Linewidth',2)
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
% plot(2.^rk,mem_app,'Linewidth',2)
% hold on
% plot(2.^rk,mem_app_TT,'--','Linewidth',2)
% xlabel('Number of particles','Fontsize',12)
% legend('Memory application TTNO binary tree','Memory application TTNO TT'...
%       ,'Fontsize',20)



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

%     %% error check of MPS
%     % build up the full operator
%     B = cell(1,d);
%     for ii=1:d
%         B{ii,ii} = Omega*sx + Delta*n_Pu;
%     end
%     count = d + 1;
%     for ii=1:d
%         for jj=1:d
%             if ii ~= jj
%                 B{count,ii} = 0.5*nu*(1/abs(ii-jj))^alpha * n_Pu;
%                 B{count,jj} = n_Pu;
%                 count = count + 1;
%             end
%         end
%     end
%     H_ex = make_operator(X_TT,B,tau_TT,2*ones(d,1));
%     H_ex = rounding(H_ex,tau_TT);
%
%     H_ex{end} = -H_ex{end};
%     D = Add_TTN(H_ex,TTNO_TT,tau_TT);
%     err = sqrt(abs(Mat0Mat0(D,D)))


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