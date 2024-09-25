function [X] = exact_Laplace(X,dt,D,d)

for ii=1:d
    mat = expm(-1i*dt*D{ii});
    X = mult_to_leaf(X,mat,ii,d);
end

end
