clear all; clc; close all;

addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\TTNO')
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\TTNO\HSS_and_no_structure_case')
addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\hm-toolbox-master')
addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\tensor_toolbox-master')

%% initializations
d = 2^9;           % number of particles
l = log(d)/log(2); % number of layers
n = 2;             % physical dimension

sx=[0,1;1,0];      % Pauli Matrix \sigma_x
n_Pu=[1,0;0,0];    %% Projector onto the excited state Pu=(sz+id)/2;
nu = 2;
Delta = 1;
Omega = 1;
alpha = 0;

[X,tau] = init_spin_all_dim_same_rank(4,1,d); % initial data - needed for construct exact operator

%% interaction matrix - single particle
A_single = cell(1,d);
for ii=1:d
    A_single{ii} = Omega*sx + Delta*n_Pu;
end
V_single = eye(d,d);

%% interaction matrix - long-range interactions
V_int = zeros(d,d);
A_int = cell(1,d);
for ii=1:d
    for jj=1:d
        if ii ~= jj
            V_int(ii,jj) = nu*(1/abs(ii-jj))^alpha;
        end
    end
    A_int{ii} = n_Pu;
end

%% HSS 
H_single = hss(V_single,'cluster',1:d);
H_single = adjust_H(H_single);

H_int = hss(V_int,'cluster',1:d);
H_int = adjust_H(H_int);

k = hssrank(H_int); 

%% construction of TTNO with help of HSS
TTNO_single = TTNO_HSS_construct_single(H_single,A_single,l,l,n*ones(1,d),1:d);
TTNO_int = TTNO_HSS_construct(H_int,A_int,l,l,n*ones(1,d),1:d);

TTNO = Add_TTN(TTNO_single,TTNO_int,tau);
TTNO = rounding(TTNO,tau);

%% exact TTNO
B = linearisation_long_range_unitary(d,sx,n_Pu,nu,Delta,Omega,alpha);
TTNO_exact_single = make_operator(X,B,tau,n*ones(1,d));

TTNO_exact_int = TTNO_no_structure(A_int,V_int,l,l,n*ones(d,1),1:d);

TTNO_exact = Add_TTN(TTNO_exact_single,TTNO_exact_int,tau);


%% error check
tmp = TTNO;
tmp{end} = -tmp{end};
E = Add_TTN(TTNO_exact,tmp,tau);
err = sqrt(abs(Mat0Mat0(E,E)))

err_scaled = sqrt(abs(Mat0Mat0(E,E)))/sqrt(abs(Mat0Mat0(TTNO_exact,TTNO_exact)))

max_rk = max_rank(TTNO)

