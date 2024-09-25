function [Y,tau] = init_Gaussian_tensor_train()
% This funciton creates the initial data of a product of Gauissians, 
% d = 2^N in the form of a binary tree
global K r d a b q N r_tau


[tau,Y] = create_tensor_train(d);

for k=1:d
    n = (-(K(k)-1)/2):1:((K(k)-1)/2);
    F_coef = zeros(K(k),1);
    for i=1:K(k)
        f = @(x) exp(-((x-q(k)).^2) / 2).*...
                 exp(-(2*pi*1i*n(i)*(x-a))/(b-a));
        F_coef(i) = (1/(b-a))*integral(f,a,b);
    end
    
    U = [F_coef (rand(K(k),r(k)-1)+1i*rand(K(k),r(k)-1))];
    Y = set_leaf_tensor_train(Y,U,k,1);
end
Y{end-1} = 1;
T = zeros(r(1),r(1));
T(1,1) = 1;
Y{end} = tensor(T,[r(1) r(1) 1]);

Y = orth_tree(Y);

end

function [tau,Y] = create_tensor_train(k)
global r_tau r d

if k==1
    tau = cell(1,2);
    Y = cell(1,4);
    Y{end-1} = eye(r(end-1),r(end-1));
    Y{end-2} = 1;
    T = zeros(r(end),r(end-1));
    T(1,1) = 1;
    Y{end} = tensor(T,[r(end) 1 r(end-1)]);
else
    tau = cell(1,2);    
    Y = cell(1,4);
    tau{1} = [];
    Y{1} = [];
    [tau{2},Y{2}] = create_tensor_train(k-1);
    if k~=d
        Y{end-1} = eye(r_tau,r_tau);
        a = d - k + 1;
        T = zeros(r(a),r(a),r(a-1));
        T(1,1,1) = 1;
        Y{end} = tensor(T,[r(a) r(a) r(a-1)]);
    end
    
end
end

function [Z] = set_leaf_tensor_train(Y,S,k,a)
% S is the matrix that shall be on k-th leaf

Z = Y;

if k==a
    Z{1} = S;
else
    Z{2} = set_leaf_tensor_train(Y{2},S,k,a+1);
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