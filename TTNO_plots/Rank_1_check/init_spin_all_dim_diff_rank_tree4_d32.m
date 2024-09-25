function [Y,tau] = init_spin_all_dim_diff_rank_tree4_d32(d)
% This funciton creates the initial data of a product state
 
% tau
tau = cell(1,4);

tau{1} = cell(1,4);
tau{1}{1} = cell(1,2);
tau{1}{2} = cell(1,2);
tau{1}{3} = cell(1,2);
tau{1}{4} = cell(1,2);

tau{2} = cell(1,4);
tau{2}{1} = cell(1,2);
tau{2}{2} = cell(1,2);
tau{2}{3} = cell(1,2);
tau{2}{4} = cell(1,2);

tau{3} = cell(1,4);
tau{3}{1} = cell(1,2);
tau{3}{2} = cell(1,2);
tau{3}{3} = cell(1,2);
tau{3}{4} = cell(1,2);

tau{4} = cell(1,4);
tau{4}{1} = cell(1,2);
tau{4}{2} = cell(1,2);
tau{4}{3} = cell(1,2);
tau{4}{4} = cell(1,2);

% Y
Y = cell(1,6);

Y{1} = cell(1,6);
Y{1}{1} = cell(1,4);
Y{1}{2} = cell(1,4);
Y{1}{3} = cell(1,4);
Y{1}{4} = cell(1,4);

Y{2} = cell(1,6);
Y{2}{1} = cell(1,4);
Y{2}{2} = cell(1,4);
Y{2}{3} = cell(1,4);
Y{2}{4} = cell(1,4);

Y{3} = cell(1,6);
Y{3}{1} = cell(1,4);
Y{3}{2} = cell(1,4);
Y{3}{3} = cell(1,4);
Y{3}{4} = cell(1,4);

Y{4} = cell(1,6);
Y{4}{1} = cell(1,4);
Y{4}{2} = cell(1,4);
Y{4}{3} = cell(1,4);
Y{4}{4} = cell(1,4);


C = zeros(2,2,4);
C(1) = 1;
C = tensor(C);
for k=1:d
    U = [1;0];
    for l=1:1
        U = [U ones(2,1)];
    end
end


Y{1}{1}{end} = C;
Y{1}{1}{end-1} = eye(4,4);
Y{1}{1}{1} = U;
Y{1}{1}{2} = U;

Y{1}{2}{end} = C;
Y{1}{2}{end-1} = eye(4,4);
Y{1}{2}{1} = U;
Y{1}{2}{2} = U;

Y{1}{3}{end} = C;
Y{1}{3}{end-1} = eye(4,4);
Y{1}{3}{1} = U;
Y{1}{3}{2} = U;

Y{1}{4}{end} = C;
Y{1}{4}{end-1} = eye(4,4);
Y{1}{4}{1} = U;
Y{1}{4}{2} = U;

C = zeros(4,4,4,4,4);
C(1) = 1;
C = tensor(C);
Y{1}{end} = C;
Y{1}{end-1} = eye(4,4);

% make the rest as Y{1}
Y{2} = Y{1};
Y{3} = Y{1};
Y{4} = Y{1};

% root tensor
C = zeros(4,4,4,4);
C(1) = 1;
Y{end} = tensor(C,[4 4 4 4 1]);
Y{end-1} = 1;

% orhtonormalization
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
