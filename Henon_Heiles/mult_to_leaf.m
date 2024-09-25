function [Y] = mult_to_leaf(Y,A,k,d)
% This function applies the matrix A to the k-th leaf of Y

if d==2
    Y{k} = A*Y{k};
    
elseif d==3 
    if k==3
        Y{2} = A*Y{2};
    else
        Y{1} = mult_to_leaf(Y{1},A,k,2);
    end
else
    N = 0;
    c = d - 2^N;
    while 2^N <= c
        N = N + 1;
        c = d - 2^N;
    end
    if c<= (2^(N-1))
        dim_left = 2^(N-1) + c;
        dim_right = 2^(N-1);
    else
        dim_left = 2^N;
        dim_right = c;
    end

    if k <= dim_left
        k_new = k;
        d_new = dim_left;
        Y{1} = mult_to_leaf(Y{1},A,k_new,d_new);
    else
        k_new = k - dim_left;
        d_new = dim_right;
        Y{2} = mult_to_leaf(Y{2},A,k_new,d_new);
    end
    
end     
end