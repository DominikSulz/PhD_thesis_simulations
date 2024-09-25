% Testscript for randomized svd

k = 10;
q = 2;
r = 5;
n = 100;
m = 200;

% initial data of rank r
U0 = orth(rand(n,r)) + 1i*orth(rand(n,r));
V0 = orth(rand(m,r)) + 1i*orth(rand(m,r));
S0 = rand(r,r) + 1i*rand(r,r);
A = U0*S0*V0';

% stage A
tic
Om = randn(m,2*k);

Y = (A*A')^q * (A*Om);

[Q,~] = qr(Y,0);

% stage B
B = Q'*A;
[U_tilde,S,V] = svd(B);
U = Q*U_tilde;
toc

tic
[U2,S2,V2] = svd(A);
toc

% test of approximation
err = norm(A - U*S*V')