function [B] = linearisation_long_range_single(d,J,sx,n,V,Delta,Omega,gamma,alpha)

B = cell(d,d);

id = [1,0;0,1];

tmp1 = gamma*(kron(J,conj(J)) - 0.5*kron(J'*J,id) - 0.5*kron(id,(J'*J).'));
tmp2 = Omega*1i*(-kron(sx,id) + kron(id,sx.'));
tmp3 = Delta*1i*(-kron(n,id) + kron(id,n.'));

for k=1:d
    B{k,k} = tmp1 + tmp2 + tmp3;
end

end