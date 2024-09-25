function [sum_rk] = sum_ranks(X)

m = length(X) - 2;
s = size(X{end});
sum_rk = sum(s(1:end-1));

for ii=1:m
    if 1==iscell(X{ii})
        sum_rk = sum_rk + sum_ranks(X{ii});
    end
end

end