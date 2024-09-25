function [TTNO] = TTNO_HSS_construct_unstructured(H,A,l,num_l,n,pos)

TTNO = cell(1,4);
I = eye(n(1),n(1));
I = I(:);

if l==1 % leaves
    tmp = A{1};
    Ui = [I tmp(:) tmp(:)];
    TTNO{1} = Ui;
    
    tmp = A{2};
    Ui = [I tmp(:) tmp(:)];
    TTNO{2} = Ui;
else
    len = length(A);
    s = len/2;
    TTNO{1} = TTNO_HSS_construct_unstructured(H.A11,A(1:s),l-1,num_l,n(1:s),pos(1:s));
    TTNO{2} = TTNO_HSS_construct_unstructured(H.A22,A((s+1):len),l-1,num_l,n((s+1):len),pos((s+1):len));
end
% connecting tensor
m1 = count_leaves(TTNO{1}); % number of leaves left
m2 = count_leaves(TTNO{2}); % number of leaves right

mat_C = zeros((m1+2)*(m2+2),m1+m2+2);

mat_C(1,1) = 1; % identity

% double interaction
% build up U and V from toolbox
[~,U,~] = full(H.A11);
[~,~,V] = full(H.A22);

Int = U * H.B12 * V'; % subdiagonal block V(\tau_1,\tau_2)

tmp = zeros(m1+2,m2+2);
if l~=1
    tmp(1,2) = 1; % h_\tau_1
    tmp(2,1) = 1; % h_\tau_2
end
tmp(3:end,3:end) = Int;
mat_C(:,2) = tmp(:); 

% single interaction
E = eye(m1+2,m2+2);
u1 = E(:,1);
for ii=1:m1
    u = E(:,ii+2);
    tmp = u*u1';
    mat_C(:,ii+2) = tmp(:); 
end

for ii=1:m2
    u = E(:,ii+2);
    tmp = u1*u';
    mat_C(:,m1+ii+2) = tmp(:);
end

% tensorization of the mat_C matrix
if l==num_l
    mat_C = mat2tens(mat_C(:,2).',[m1+2 m2+2 1],3,1:2);
    TTNO{end-1} = 1;
    TTNO{end} = tensor(mat_C,[m1+2 m2+2 1]);
else
    mat_C = mat2tens(mat_C.',[m1+2 m2+2 m1+m2+2],3,1:2);
    TTNO{end-1} = eye(m1+m2+2,m1+m2+2);
    TTNO{end} = tensor(mat_C,[m1+2 m2+2 m1+m2+2]);
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