function [B] = linearisation_Ising(sx,sz,d)
% This function saves the matrices, which are needed to represent the
% discrete Hamiltonian, in an array B


B = cell(d,d);
for ii=1:d
    B{ii,ii} = sx;
end

count = d+1;
for ii=1:(d-1)
    B{count,ii} = sz;
    B{count,ii+1} = sz;
    count = count + 1;
end

end
