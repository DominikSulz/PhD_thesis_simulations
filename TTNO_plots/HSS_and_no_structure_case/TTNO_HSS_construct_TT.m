function [TTNO] = TTNO_HSS_construct_TT(H,A,l,num_l,n,pos)

TTNO = cell(1,4);
I = eye(n(1),n(1));
I = I(:);

if length(pos) == 2 % leafs 
    tmp = A{1};
    Ui = [I tmp(:) tmp(:)];
    TTNO{1} = Ui;
    
    tmp = A{2};
    Ui = [I tmp(:) tmp(:)];
    TTNO{2} = Ui;
else
    tmp = A{1};
    Ui = [I tmp(:) tmp(:)];
    TTNO{1} = Ui;

    len = length(A);
    TTNO{2} = TTNO_HSS_construct_TT(H.A22,A(2:len),l-1,num_l,n(2:len),pos(2:len));
end

%% connecting tensor
[~,U,~] = full(H.A11);
[~,~,V] = full(H.A22);
[~,t1] = size(U);
[~,t2] = size(V);
ss = max(t1,t2);

m1 = count_leaves(TTNO{1});
m1 = min([m1,t1]);
m2 = count_leaves(TTNO{2});
m2 = min([m2,t2]);

[~,U,~] = full(H);
[~,ss2] = size(U);
mat_C = zeros((m1+2)*(m2+2),min([m1+m2,ss2])+2); % \mat_{1:2}(C)

%% identity
mat_C(1,1) = 1; 

%% single interaction 
single = [];
E = eye(m1+2,m2+2);
u1 = zeros(m2+2,1);
u1(1) = 1;
for ii=1:m1          % u_i*u_1^T, where i=3,...,2+k
%     u = E(ii+2,:)';
    u = zeros(m1+2,1);
    u(ii+2) = 1;
    tmp = u*u1';
    single(:,ii) = tmp(:); 
end

u1 = zeros(m1+2,1);
u1(1) = 1;
for ii=1:m2         % u_1*u_i^T, where i=3,...,2+k
    u = zeros(m2+2,1);
    u(ii+2) = 1;
%     u = E(:,ii+2);
    tmp = u1*u';
    single(:,m1+ii) = tmp(:);
end

if l ~= num_l 
    R = [H.Rl;H.Rr];
        
    if isempty(single) == 0
        mat_C(:,3:end) = single*R; 
    end
end

%% double interaction
tmp = zeros(m1+2,m2+2);
if l~=1 % give h_\tau_1 and h_\tau_2 one level up
    tmp(1,2) = 1; 
    % tmp(2,1) = 1; % deleted, as we dont give h_\tau_1 = A_i to the upper
    %level twice
end

tmp(3:end,3:end) = H.B12; % subdiagonal block V(\tau_1,\tau_2)
mat_C(:,2) = tmp(:); 

if l==num_l
    mat_C = mat2tens(mat_C(:,2).',[m1+2 m2+2 1],3,1:2);
    TTNO{end-1} = 1;
    TTNO{end} = tensor(mat_C,[m1+2 m2+2 1]);
else
    mat_C = mat2tens(mat_C.',[m1+2 m2+2 min([m1+m2,ss2])+2],3,1:2);
    TTNO{end-1} = eye(min(m1+m2,ss2)+2,min([m1+m2,ss2])+2);
    TTNO{end} = tensor(mat_C,[m1+2 m2+2 min([m1+m2,ss2])+2]);
end


end


function [erg] = count_leaves(X)

if iscell(X) == 0
    erg = 1;
else
    m = length(X) - 2;
    arr = zeros(m+1,1);
    
    for ii=1:m
        if 1==iscell(X{ii})
            tmp = count_leaves(X{ii});
            arr(ii) = tmp(end);
        else
            arr(ii) = 1;
        end
    end
    arr(end) = sum(arr(1:end-1));
    erg = arr(end);
end
end


%%% test for single interaction
%     %%%% test
% %     [u,s,v] = svd(R,'econ');
%     [~,ss] = size(R);
%     [rr,~] = size(mat_C);
%     if ss>k
%         mat_C(:,(end+1):(end+ss-k)) = zeros(rr,ss-k);
%         k = ss;
%     end
%     % if I do that, then an error comes where sometimes U_\tau_i has less
%     % than n_i rows - does also not fit theory
%     
%     %%%%


%%%% old code


% % single interaction (first idea)
% [~,U1,V1] = full(H.A11);
% [~,U2,V2] = full(H.A22);
% % U = zeros(size(U1)+size(V2));
% [s1,s2] = size(U1);
% % U(1:s1,1:s2) = U1;
% [t1,t2] = size(V2);
% % U((s1+1):(s1+t1),(s2+1):(s2+t2)) = V2;
% % R = [H.Rl;H.Rr];
% % UR = U*R;
% % [~,s_R] = size(UR);
% 
% [~,U,V] = full(H);
% n_l1 = count_leaves(TTNO{1});
% n_l2 = count_leaves(TTNO{2});
% E = eye(n_l1+2,n_l2+2);
% u1 = E(:,1);
% for ii=1:n_l1
%     u = E(:,ii+2);
%     tmp = u*u1';
%     tmp(3:end,3:end) = U'*tmp(3:end,3:end);
%     K(:,ii) = tmp(:); 
% end
% 
% for ii=1:n_l2
%     u = E(:,ii+2);
%     tmp = u1*u';
%     tmp(3:end,3:end) = U'*tmp(3:end,3:end);
%     K(:,m1+ii) = tmp(:);
% end
% mat_C(:,3:end) = K;

% % single interaction (idea 2)
% if l~=num_l
%     
%     if (m1+m2) <= k
%         E = eye(m1+2,m2+2);
%         u1 = E(:,1);
%         for ii=1:m1
%             u = E(:,ii+2);
%             tmp = u*u1';
%             mat_C(:,ii+2) = tmp(:);
%         end
%         
%         for ii=1:m2
%             u = E(:,ii+2);
%             tmp = u1*u';
%             mat_C(:,m1+ii+2) = tmp(:);
%         end
%     else
%         R = [H.Rl;H.Rr];
%         for ii=1:min(m1+m2,k)
%             tmp = R(:,ii);
%             mat_C(3:end,2+ii) = tmp(:);
%         end
%     end
% end