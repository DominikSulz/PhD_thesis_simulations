function [TTNO] = TTNO_HSS_construct_single_TT(H,A,l,num_l,n,pos)

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
    TTNO{2} = TTNO_HSS_construct_single_TT(H.A22,A(2:len),l-1,num_l,n(2:len),pos(2:len));
end

%% connecting tensor

mat_C = zeros(9,3);
mat_C(1,1) = 1; % identity
mat_C(2,2) = 1; % left part
mat_C(4,2) = 1; % right part % davor (3,2)

if l==num_l
    mat_C = mat2tens(mat_C(:,2).',[3 3 1],3,1:2);
    TTNO{end-1} = 1;
    TTNO{end} = tensor(mat_C,[3 3 1]);
else
    mat_C = mat2tens(mat_C.',[3 3 3],3,1:2);
    TTNO{end-1} = eye(2,2);
    TTNO{end} = tensor(mat_C,[3 3 3]);
end


end