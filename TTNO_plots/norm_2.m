function[nn] = norm_2(H,samp,d)
% approximates the 2 norm of an operator.

r = 3;

nn_save = zeros(1,samp);
for jj=1:samp
    X = rand_tree(H,r);
    % normalize
    X{end} = X{end}/sqrt(abs(Mat0Mat0(X,X)));
    % apply operator
    tmp = apply_operator_nonglobal(X,H,d);
    nn_save(jj) = sqrt(abs(Mat0Mat0(tmp,tmp)));
end

nn = max(nn_save);

end



function[X] = rand_tree(H,r)

X = H;
m = length(X) - 2;
n = size(H{end});

if n(end) == 1
    tmp = randn(r*ones(1,m));
    X{end} = tensor(tmp,[r*ones(1,m) 1]);
else
    tmp = randn(r*ones(1,m+1));
    X{end} = tensor(tmp,r*ones(1,m+1));
end


for ii=1:m
    if iscell(X{ii}) == 1
        X{ii} = rand_tree(H{ii},r);
    else
        [n,~] = size(H{ii});
        X{ii} = randn(sqrt(n),r);
    end
end


end