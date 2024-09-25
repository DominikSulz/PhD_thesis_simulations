function [Y,tau] = init_spin_all_dim_diff_rank_tree4(d)
% This funciton creates the initial data of a product state
 
tau = cell(1,4);
tau{1} = cell(1,4);
tau{2} = cell(1,4);
tau{3} = cell(1,4);
tau{4} = cell(1,4);

Y = cell(1,6);
Y{1} = cell(1,6);
Y{2} = cell(1,6);
Y{3} = cell(1,6);
Y{4} = cell(1,6);
C = zeros(2,2,2,2,4);
C(1) = 1;
C = tensor(C);
Y{1}{end} = C;
Y{1}{end-1} = eye(4,4);
Y{2}{end} = C;
Y{2}{end-1} = eye(4,4);
Y{3}{end} = C;
Y{3}{end-1} = eye(4,4);
Y{4}{end} = C;
Y{4}{end-1} = eye(4,4);
C = zeros(4,4,4,4);
C(1) = 1;
C = tensor(C,[4 4 4 4 1]);
Y{end} = C;
Y{end-1} = 1;

for k=1:d
    U = [1;0];
    for l=1:1
        U = [U ones(2,1)];
    end
end

for ii=1:d
    if ii<=4
        Y{1}{ii} = U;
    elseif 4 < ii && ii<= 8
        Y{2}{ii-4} = U;
    elseif 8 < ii && ii<= 12
        Y{3}{ii-8} = U;
    elseif 12 < ii && ii<= 16
        Y{4}{ii-12} = U;
    end
end

Y = orth_tree(Y);

end

function [X] = orth_tree(T)
% This function orthonormalizes the core of T

X = T;
m = length(T) - 2;

for i=1:m
    if 1==iscell(T{i})
        dum2 = orth_tree(T{i});
        s = size(dum2{end});
        ss = length(s);
        Mat = double(tenmat(dum2{end},ss,1:(ss-1)));
        [Q,R] = qr(Mat.',0);
        dum2{end} = tensor(mat2tens(Q.',s,ss,1:(ss-1)));
        X{i} = dum2;
        X{end} = ttm(X{end},R,i);
    else
        [Q,R] = qr(T{i},0);
        X{i} = Q;
        X{end} = ttm(X{end},R,i);
    end
end

end
