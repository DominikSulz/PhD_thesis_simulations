function [Y,tau,F_coef] = init_Gaussian_tucker()
% This funciton creates the initial data of a product of Gauissians
global K r d a b q

tau = cell(1,d);

Y = cell(1,d+2);

Y{end-1} = 1;
T = zeros([r 1]);
l = ones(length(r)+1,1);
T(l) = 1;
Y{end} = tensor(T,[r 1]);

for k=1:d
    n = (-(K(k)-1)/2):1:((K(k)-1)/2);
    F_coef = zeros(K(k),1);
    for i=1:K(k)

        f = @(x) exp(-((x-q(k)).^2) / 2).*...
                 exp(-(2*pi*1i*n(i)*(x-a))/(b-a));
        F_coef(i) = (1/(b-a))*integral(f,a,b);
    end
    
    U = [F_coef (rand(K(k),r(k)-1)+1i*rand(K(k),r(k)-1))];
    [Q,R] = qr(U,0);
    Y{k} = Q;
    Y{end} = ttm(Y{end},R,k);
    
end

end