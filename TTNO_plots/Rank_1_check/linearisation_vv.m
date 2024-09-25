function [B] = linearisation_vv(A,d,v)
% linearisation of H = \sum_{ii<jj}^d V(i,j)A^(i)A^(j),
% where V(i,j) = v(i)*v(j), for some function/vector v.

B = cell(1,d);
count = 1;

for ii=1:d
    for jj=1:d
        if ii<jj
            B{count,ii} = v(ii)*A;
            B{count,jj} = v(jj)*A;
            count = count + 1;
        end
    end
end

end