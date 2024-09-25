function [Y,tau] = init_spin_all_dim_diff_rank_TT(d)
% This funciton creates the initial data of a product state as 

[tau,Y] = create_TT(d);

for k=1:d
    U = [1;0];
    for l=1:1
        U = [U ones(2,1)];
    end
    Y = set_leaf_TT(Y,U,k,d);
end
Y{end-1} = 1;

C = zeros(2,2);
C(1) = 1;
Y{end} = tensor(C,[2 2 1]);

if d>2
    C = zeros(2,4,2);
    C(1) = 1;
    Y{2}{end} = tensor(C,[2 4 2]);
end

Y = orth_tree(Y);

end

function [tau,Y] = create_TT(d)

if d==2
    tau = cell(1,3);
    Y = cell(1,4);
else
    tau = cell(1,2);
    Y = cell(1,4);
    [tau{2},Y{2}] = create_TT(d-1);
end


end


function [Z] = set_leaf_TT(Y,S,k,d)
% S is the matrix that shall be on k-th leaf

Z = Y;
if d==2 && k<=2 
    Z{k} = S;
    C = zeros(2,2,4);
    C(1) = 1;
    Z{end} = tensor(C,[2 2 4]);
    Z{end-1} = eye(4,4);
elseif k==1 && d>=3
    C = zeros(2,4,4);
    C(1) = 1;
    Z{end} = tensor(C,[2 4 4]);
    Z{end-1} = eye(4,4);
    Z{1} = S;
else
    Z{2} = set_leaf_TT(Y{2},S,k-1,d-1);
end

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


%% Old
% function [Z] = set_leaf_TT(Y,S,k,d)
% % S is the matrix that shall be on k-th leaf
% 
% Z = Y;
% if d==1 
%     Z{1} = S;
% %     Z{2} = 1;
%     C = zeros(2,2,4);
%     C(1) = 1;
%     Z{end} = tensor(C,[2 2 4]);
%     Z{end-1} = eye(2,2);
% elseif k==1
%     Z{k} = S;
% else
%     if d==2
%         C = zeros(2,2,4);
%         C(1) = 1;
%         Z{end} = tensor(C,[2 2 4]);
%         Z{end-1} = eye(4,4);
%         Z{2} = set_leaf_TT(Y{2},S,k-1,d-1);
%     else
%         C = zeros(2,4,4);
%         C(1) = 1;
%         Z{end} = tensor(C,[2 4 4]);
%         Z{end-1} = eye(4,4);
%         Z{2} = set_leaf_TT(Y{2},S,k-1,d-1);
%     end
% end
% 
% end