function [Y] = rounding2(A,tau)
Y = A;
m = length(A) - 2;
for i=1:m
    if 1 == iscell(A{i}) % A{i} is a TTN
        dum = rounding2(A{i},tau{i});
        [Y{i},Ri] = orth_core(dum);        
        Y{end} = ttm(Y{end},Ri,i); % orthonormalize the core
    else % if we are at a basis matrix
        [Qi,S,V] = svd(A{i},"econ"); % davor ohne econ
        Ri = S*V';             
%         [Qi,Ri] = qr(A{i},0);     % when truncating its better to use
%         svd, as we want the maximum eigenvalues
        
        [s1,s2] = size(S);
        ss = min(s1,s2);
        tol_S = eps;%*norm(S,'Fro'); %Änderung: davor gar nicht 
        
        S_diag = diag(S);
        rk = [];
        for j=1:ss
            tmp = sqrt(sum(S_diag(j:ss).^2));
            if tmp < tol_S % zuvor tol_S/m
                rk = j-1;
                break
            end
        end
        if 1==isempty(rk)
            rk = ss;
        end
        rk = max(rk,2);
        rk = min(rk,100);
        r_tilde = rk;
        
        if r_tilde ~= 1
            Qi = Qi(:,1:r_tilde);
            Ri = Ri(1:r_tilde,:);
        end
        Y{i} = Qi;
        Y{end} = ttm(Y{end},Ri,i); 
    end
end
s = size(Y{end});
Y{end-1} = eye(s(end),s(end));

end

function [X,R] = orth_core(T)
% This function orthonormalizes the core of T

m = length(T) - 2;
Mat = double(tenmat(T{end},m+1,1:m));
[Mat,S,V] = svd(Mat.','econ'); % davor ohne econ

R = S*V';
s = size(T{end});

[s1,s2] = size(S);
ss = min(s1,s2);
tol_S = eps;%*norm(S,'Fro'); %Änderung: davor gar nicht

S_diag = diag(S);
rk = [];
for j=1:ss
    tmp = sqrt(sum(S_diag(j:ss).^2));
    if tmp < tol_S % zuvor tol_S/m
        rk = j-1;
        break
    end
end
if 1==isempty(rk)
    rk = ss;
end
rk = max(rk,2);
rk = min(rk,100);
rC = rk;
% 
% rC = rank(R);

if rC ~= 1
    Mat = Mat(:,1:rC);
    R = R(1:rC,:);
    s(end) = rC;
end
if size(Mat.',1) ~= prod(s(m+1)) || ...
   size(Mat.',2) ~= prod(s(1:m))
%     1
end
X = T;
X{end} = tensor(mat2tens(Mat.',s,m+1,1:m),s);
X{end-1} = eye(rC,rC);
end