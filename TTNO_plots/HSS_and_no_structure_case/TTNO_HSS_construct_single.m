function [TTNO] = TTNO_HSS_construct_single(H,A,l,num_l,n,pos)

TTNO = cell(1,4);
I = eye(n(1),n(1));
I = I(:);

if l==1 % leafs
    tmp = A{1};
    Ui = [I tmp(:)];
    TTNO{1} = Ui;
    
    tmp = A{2};
    Ui = [I tmp(:)];
    TTNO{2} = Ui;
else
    % recursion - build up the subtrees \tau_1 and \tau_2
    len = length(A);
    s = len/2;
    TTNO{1} = TTNO_HSS_construct_single(H.A11,A(1:s),l-1,num_l,n(1:s),pos(1:s));
    TTNO{2} = TTNO_HSS_construct_single(H.A22,A((s+1):len),l-1,num_l,n((s+1):len),pos((s+1):len));
end

%% connecting tensor

mat_C = zeros(4,2);
mat_C(1,1) = 1; % identity
mat_C(2,2) = 1; % left part
mat_C(3,2) = 1; % right part

if l==num_l
    mat_C = mat2tens(mat_C(:,2).',[2 2 1],3,1:2);
    TTNO{end-1} = 1;
    TTNO{end} = tensor(mat_C,[2 2 1]);
else
    mat_C = mat2tens(mat_C.',[2 2 2],3,1:2);
    TTNO{end-1} = eye(2,2);
    TTNO{end} = tensor(mat_C,[2 2 2]);
end


end