function [sum_rk,N] = prod_ranks(X)

m = length(X) - 2;
sum_rk = prod(size(X{end}));
N = 1;

for ii=1:m
    if iscell(X{ii}) == 1
        [rk_i,N_i] = prod_ranks(X{ii});
        sum_rk = sum_rk + rk_i;
        N = N + N_i;
    end
end

end