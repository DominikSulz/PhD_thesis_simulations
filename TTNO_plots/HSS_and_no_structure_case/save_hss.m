function [TTNO] = TTNO_HSS_construct(H,A,l,num_l,n,pos,k)
%% for now just copied - has to be adjusted
TTNO = cell(1,4);
I = eye(n(1),n(1));
I = I(:);

if l==1 % leaf 
    tmp = A{1};
    Ui = [I tmp(:) tmp(:)];
    TTNO{1} = Ui;
    
    tmp = A{2};
    Ui = [I tmp(:) tmp(:)];
    TTNO{2} = Ui;
else
    len = length(A);
    s = len/2;
    TTNO{1} = TTNO_HSS_construct(H.A11,A(1:s),l-1,num_l,n(1:s),pos(1:s),k);
    TTNO{2} = TTNO_HSS_construct(H.A22,A((s+1):len),l-1,num_l,n((s+1):len),pos((s+1):len),k);
end
% connecting tensor
kk = hssrank(H);
m1 = count_leaves(TTNO{1});
m1 = min([m1,k]);
m2 = count_leaves(TTNO{2});
m2 = min([m2,k]);

mat_C = zeros((m1+2)*(m2+2),min([m1+m2,k])+2); % old mat_C = zeros((m1+2)*(m2+2),min(m1+m2,k)+2);

mat_C(1,1) = 1; % identity

% single interaction 
E = eye(m1+2,m2+2);
u1 = E(:,1);
for ii=1:m1
    u = E(:,ii+2);
    tmp = u*u1';
    single(:,ii) = tmp(:); 
end

for ii=1:m2
    u = E(:,ii+2);
    tmp = u1*u';
    single(:,m1+ii) = tmp(:);
end
if l ~= num_l
    mat_C(:,3:end) = single*[H.Rl;H.Rr];
end

% double interaction
tmp = zeros(m1+2,m2+2);
if l==1
    
else
    tmp(1,2) = 1;
    tmp(2,1) = 1;
end

tmp(3:end,3:end) = H.B12; % subdiagonal block V(\tau_1,\tau_2)
mat_C(:,2) = tmp(:); 

if l==num_l
    mat_C = mat2tens(mat_C(:,2).',[m1+2 m2+2 1],3,1:2);
    TTNO{end-1} = 1;
    TTNO{end} = tensor(mat_C,[m1+2 m2+2 1]);
else
%     [~,ss] = size(KU);
    mat_C = mat2tens(mat_C.',[m1+2 m2+2 min([m1+m2,k])+2],3,1:2);
    TTNO{end-1} = eye(min(m1+m2,k)+2,min([m1+m2,k])+2);
    TTNO{end} = tensor(mat_C,[m1+2 m2+2 min([m1+m2,k])+2]);
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