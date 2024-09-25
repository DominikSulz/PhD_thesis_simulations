function [C] = core_tensor(X)

m = length(X) - 2;

if iscell(X{1}) == 0
    mat_C = zeros(2^m,3);
else
    mat_C = zeros(3^m,3);
end

% identity
mat_C(1,1) = 1;

% vec(A^(i))
if iscell(X{1}) == 0
    for ii=1:m
        mat_C(2^(m-ii) + 1,2) = 1;
    end
else
    for ii=1:m
        mat_C(3^(m-ii) + 1,2) = 1;
    end
end

% vec(A^(i))vec(A^(j))
if iscell(X{1}) == 0
    for ii=1:m
        for jj=1:m
            if ii < jj
                mat_C(2^(m-ii) + 2^(m-jj) + 1,3) = 1;
            end
        end
    end
else
    for ii=1:m
        for jj=1:m
            mat_C(3^(m-ii) + 3^(m-jj) + 1,3) = 1;
        end
    end
end

% set core
if iscell(X{1}) == 0
    s = 2*ones(1,m);
else
    s = 3*ones(1,m);
end
s = [s 3];

tmp = size(X{end});
if tmp(end) == 1
   mat_C = mat_C(:,3); % if we are at the core tensor
   s(end) = 1;
end

C = mat2tens(mat_C.',s,m+1,1:m);
C = tensor(C,s);

end

