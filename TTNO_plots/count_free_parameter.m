function[num] = count_free_parameter(X)

m = length(X) - 2;

% count free parameters in the connecting tensor
tmp = double(X{end});
num = numel(tmp);

for ii=1:m
    if 1 == iscell(X{ii})
        num = num + count_free_parameter(X{ii});
    else
        num = num + numel(X{ii});
    end
end