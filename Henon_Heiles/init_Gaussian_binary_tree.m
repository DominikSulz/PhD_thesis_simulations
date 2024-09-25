function [Y,tau] = init_Gaussian_binary_tree()
% This funciton creates the initial data of a product of Gauissians, 
% d = 2^N in the form of a binary tree
global K r d a b q N r_tau

tau = cell(1,2);
Y = cell(1,4);
for k=1:2
    [tau{k},Y{k}] = create_binary_tree(N);
end


% Y{end-1} = 1;
% T = zeros([r 1]);
% l = ones(length(r)+1,1);
% T(l) = 1;
% Y{end} = tensor(T,[r 1]);

for k=1:d
    n = (-(K(k)-1)/2):1:((K(k)-1)/2);
    F_coef = zeros(K(k),1);
    for i=1:K(k)
        f = @(x) exp(-((x-q(k)).^2) / 2).*...
                 exp(-(2*pi*1i*n(i)*(x-a))/(b-a));
        F_coef(i) = (1/(b-a))*integral(f,a,b);
    end
    
    U = [F_coef (rand(K(k),r(k)-1)+1i*rand(K(k),r(k)-1))];
    Y = set_leaf(Y,U,k,d);
    
%     [Q,R] = qr(U,0);
%     Y{k} = Q;
%     Y{end} = ttm(Y{end},R,k);
end
Y{end-1} = 1;
T = zeros(r_tau,r_tau);
T(1,1) = 1;
Y{end} = tensor(T,[r_tau r_tau 1]);

Y = orth_tree(Y);

end

function [tau,Y] = create_binary_tree(k)

if k==1
    tau = [];
    Y = [];
else
    tau = cell(1,2);
    Y = cell(1,4);
    for i=1:2
        [tau{i},Y{i}]  = create_binary_tree(k-1);
    end
    
end
end

function [Z] = set_leaf(Y,S,k,d)
% S is the matrix that shall be on k-th leaf
global r_tau r
Z = Y;
if d==2
    Z{k} = S;
    if mod(k,2) == 0 % even k
        Z{end-1} = eye(r_tau,r_tau);
        M = zeros(r(k-1),r(k),r_tau);
        M(1,1,1) = 1;
        Z{end} = tensor(M);
    else % odd k
        Z{end-1} = eye(r_tau,r_tau);
        M = zeros(r(k),r(k+1),r_tau);
        M(1,1,1) = 1;
        Z{end} = tensor(M);
    end

else
    d_new = d/2;
    if k <= (d/2)
        k_new = k;
        Z{1} = set_leaf(Y{1},S,k_new,d_new);
    else
        k_new = k - d_new;
        Z{2} = set_leaf(Y{2},S,k_new,d_new);
    end
    M = eye(r_tau,r_tau);
    Z{end-1} = M;
    M = zeros(r_tau,r_tau,r_tau);
    M(1,1,1) = 1;
    Z{end} = tensor(M);
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