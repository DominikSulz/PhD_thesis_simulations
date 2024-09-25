function [TTNO] = rank1_TTNO_vv(X,v,A,I,d)

TTNO = X;
% core tensors
TTNO = set_TTNO_core(TTNO);

% leaves
for ii=1:d
    tmp = v(ii)*A;
    U_i = [I(:),tmp(:)];
    
    TTNO = set_operator(TTNO,U_i,ii);
end


end

function [X] = set_TTNO_core(X)

m = length(X) - 2;
X{end} = core_tensor(X);

for ii=1:m
    if iscell(X{ii}) == 1
        X{ii} = set_TTNO_core(X{ii});
    end
end
    

end

function [A] = set_operator(Y,Q,k)
% S is the matrix that shall be on k-th leaf

A = Y;
arr = count_leaves(Y);
m = length(Y) - 2;
arr2 = zeros(m,1);
arr2(1) = arr(1);
for ii=2:m
    arr2(ii) = arr2(ii-1) + arr(ii);
end

tmp = 0;
for ii=1:m
    if (k <= arr2(ii)) && (tmp < k)
        if 1 == iscell(Y{ii})
            A{ii} = set_operator(Y{ii},Q,k - tmp);
        else
            A{ii} = Q;
        end
    end
    tmp = tmp + arr(ii);
end

end

function [arr] = count_leaves(X)

m = length(X) - 2;
arr = zeros(m+1,1);

for ii=1:m
    if 1==iscell(X{ii})
        tmp = count_leaves(X{ii});
        arr(ii) = tmp(end);
    else
        arr(ii) = 1;
    end
end
arr(end) = sum(arr(1:end-1));

end