clear all; clc; close all;

addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\TTNO')
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\TTNO\HSS_and_no_structure_case')
addpath('C:\Users\Dominik\Documents\MATLAB\Low rank approximations\TTNO\example_unitary_dynamics_spin')
addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\hm-toolbox-master')
addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\tensor_toolbox-master')

d = 2^7;

V = zeros(d,d);

f = @(x) 1/sin(x);
for ii=1:d
    for jj=1:d
        V(ii,jj) = f(abs(ii-jj));
    end
end

V = 0.5*(V+V');

H = hss(V,'cluster',1:d);
k = hssrank(H)