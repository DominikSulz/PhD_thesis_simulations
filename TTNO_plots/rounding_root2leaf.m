function [Y] = rounding_root2leaf(A,tau)
Y = A;
m = length(A) - 2;

for ii=1:m
    vv = 1:(m+1);
    vv = vv(vv~=ii);
    Mat = double(tenmat(Y{end},ii,vv)).';
    [U,S,V] = svd(Mat,"econ");
    R = S*V';
    rk = rank(R);
    s = size(Y{end});
    if rk ~= 1
        U = U(:,1:rk);
        R = R(1:rk,:);
        s(ii) = rk;
    end
    Y{end} = tensor(mat2tens(U.',s,ii,vv),s);
    Y{end-1} = eye(s(end),s(end));
    
    if iscell(Y{ii}) == 1
        mm = length(Y{ii}) - 2;
        Y{ii}{end} = ttm(Y{ii}{end},R,mm+1);
        Y{ii} = rounding_root2leaf(Y{ii},tau{ii});
    else
        Y{ii} = Y{ii}*R.';
    end
end
end
