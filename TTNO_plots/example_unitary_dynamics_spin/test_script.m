clear all; clc;

addpath('C:\Users\Dominik\Documents\MATLAB\Matlab toolboxes\hm-toolbox-master')

d = 64;
B = randn(d);

%% of HSS 

cluster = hss(); 
cluster = cluster_rec_HSS(cluster,d);
cluster.topnode = 1;

H = hss(B, 'cluster', cluster);

norm(full(H) - B, 'fro')

spy(H)

function[cluster] = cluster_rec_HSS(cluster,d)

cluster.leafnode = 0;
cluster.topnode = 0;
cluster.ml = 1; cluster.nl = 1;
cluster.mr = d-1; cluster.nr = d-1;

cluster.A11 = hss(); 
cluster.A11.D = 0;
cluster.A11.leafnode = 1;
cluster.A11.topnode = 0;
cluster.A22.topnode = 0;

if d==2
    cluster.A22 = hss(); 
    cluster.A22.D = 0;
    cluster.A22.leafnode = 1;
else
    d = d-1;
    cluster.A22 = hss(); cluster.A22 = cluster_rec_HSS(cluster.A22,d);
    cluster.A22.D = zeros(d);
end

end


% %% for hodlr matrices
% % first level
% cluster = hodlr(); cluster.sz = [d d];
% cluster = cluster_rec(cluster,d);
% 
% % % first test
% % cluster.A11 = hodlr(); cluster.A11.sz = [1 1];
% % cluster.A22 = hodlr(); cluster.A22.sz = [d-1 d-1];
% % % next level
% % cluster.A22.sz = [d-1 d-1];
% % cluster.A22.A11 = hodlr(); cluster.A22.A11.sz = [1,1]; 
% % cluster.A22.A22 = hodlr(); cluster.A22.A22.sz = [d-2 d-2];
% 
% 
% H = hodlr(B, 'cluster', cluster);
% 
% norm(full(H) - B, 'fro')
% 
% spy(H)

% function[cluster] = cluster_rec(cluster,d)
% 
% cluster.A11 = hodlr(); cluster.A11.sz = [1 1];
% if d==2
%     cluster.A22 = hodlr(); cluster.A22.sz = [1 1];
% else
%     d = d-1;
%     cluster.A22 = hodlr(); cluster.A22 = cluster_rec(cluster.A22,d);
%     cluster.A22.sz = [d d];
% end
% 
% end