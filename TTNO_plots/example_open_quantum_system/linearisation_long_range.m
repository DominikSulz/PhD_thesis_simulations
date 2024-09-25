function [B] = linearisation_long_range(d,J,sx,n,V,Delta,Omega,gamma,alpha)

B = cell(d,d);

id = [1,0;0,1];

tmp1 = gamma*(kron(J,conj(J)) - 0.5*kron(J'*J,id) - 0.5*kron(id,(J'*J).'));
tmp2 = Omega*(-1i*kron(sx,id) + 1i*kron(id,sx.'));
tmp3 = Delta*(-1i*kron(n,id) + 1i*kron(id,n.'));

for k=1:d
    B{k,k} = tmp1 + tmp2 + tmp3;
end

count = d + 1;
% long range
for ii=1:d
    for jj=1:d
        if ii==jj
            
        else
            B{count,ii} = -1i*(V/2) * 1/abs(ii-jj)^alpha * kron(n,id);
            B{count,jj} = kron(n,id);
            count = count + 1;
        end
    end
end

for ii=1:d
    for jj=1:d
        if ii==jj
            
        else
            B{count,ii} = 1i*(V/2) * 1/abs(ii-jj)^alpha * kron(id,n.');
            B{count,jj} = kron(id,n.');
            count = count + 1;
        end
    end
end
end