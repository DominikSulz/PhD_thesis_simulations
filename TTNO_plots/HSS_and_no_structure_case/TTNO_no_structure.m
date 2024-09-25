function [TTNO] = TTNO_no_structure(A,V,l,num_l,n,pos)
% This code works for long-range operators of the form
% \sum_i<j A^(i)A^(j)

I = eye(n(1),n(1));
I = I(:);
TTNO = cell(1,4);

if l==1
    %% leafs of the TTNO
    tmp = A{1};
    U = [I tmp(:)];
    TTNO{1} = U;
    
    tmp = A{2};
    U = [I tmp(:)];
    TTNO{2} = U;
    
    if num_l ~= 1
        %% connecting tensor at the leafs
        mat_C = zeros(4,4);
        mat_C(1,1) = 1;  % identity
        mat_C(2,2) = 1;  % A_1
        mat_C(3,3) = 1;  % A_2
        mat_C(4,4) = V(pos(1),pos(2));  % A_1 \otimes A_2
        
        TTNO{3} = eye(4,4);
        TTNO{4} = mat2tens(mat_C.',[2 2 4],1:2,3); % mat2tens(mat_C.',[2 2 4],1:2,3); 
        TTNO{4} = tensor(TTNO{4},[2 2 4]);
    else % matrix case
        TTNO{3} = 1;
        C = zeros(2,2);
        C(2,2) = V(1,2);
        TTNO{4} = tensor(C,[2 2 1]);
    end
    
elseif (l>1) && (l<num_l)
    % core tensors in between
    len = length(A);
    s = len/2;
    TTNO{1} = TTNO_no_structure(A(1:s),V,l-1,num_l,n(1:s),pos(1:s));
    TTNO{2} = TTNO_no_structure(A((s+1):len),V,l-1,num_l,n((s+1):len),pos((s+1):len));
    
    rk = (2^(l-1)+2)^2;
    mat_C = zeros(rk,2^l+2);
    
    mat_C(1,1) = 1;    % identity
    
    count = 2; % start at colum 2, as first is identities
    for jj=1:(2^(l-1))
         mat_C(jj+1,count) = 1;    % A_j (from the left subtree)
         count = count + 1;
    end
    
    for jj=1:(2^(l-1))
        tmp = jj*(2^(l-1) + 2) + 1;
        mat_C(tmp,count) = 1;    % A_j (from the right subtree)
        count = count + 1;
    end
    
    
    mat_C(2^(l-1)+2,2^l+2) = 1;     % A_i \otimes A_j which are only at the right subtree
    tmp = (2^(l-1) + 2)^2 - (2^(l-1) + 2) + 1;
    mat_C(tmp,2^l+2) = 1;           % A_i \otimes A_j which are only at the left subtree
    
    for ii=1:(2^(l-1))
        for jj=1:(2^(l-1))
            mat_C(jj*(2^(l-1)+2) + ii + 1,2^l+2) = V(pos(ii),pos(2^(l-1)+jj));  % A_i \otimes A_j, A_i from left and A_j from right subtree
        end
    end
   
    TTNO{3} = eye(2^l+2,2^l+2);
    TTNO{4} = mat2tens(mat_C.',[2^(l-1)+2 2^(l-1)+2 2^l+2],3,1:2);
    TTNO{4} = tensor(TTNO{4},[2^(l-1)+2 2^(l-1)+2 2^l+2]);
else
    % root tensor
    len = length(A);
    s = len/2;
    TTNO{1} = TTNO_no_structure(A(1:s),V,l-1,num_l,n(1:s),pos(1:s));
    TTNO{2} = TTNO_no_structure(A((s+1):len),V,l-1,num_l,n((s+1):len),pos((s+1):len));
    
    TTNO{end-1} = 1;
    mat_C = zeros((2^(l-1)+2)^2,1);
    
    for ii=1:(2^(l-1))
        for jj=1:(2^(l-1))
            mat_C(jj*(2^(l-1)+2) + ii +1) = V(pos(ii),pos(2^(l-1)+jj));  % A_i \otimes A_j, A_i from left and A_j from right subtree
        end
    end
    
    mat_C(2^(l-1)+2) = 1;                           % \sum Ai*Aj from right subtree
    mat_C((2^(l-1)+2)^2 - (2^(l-1) + 2) + 1) = 1;   % \sum Ai*Aj from left subtree
%     mat_C = reshape(mat_C,[2^(l-1)+2 2^(l-1)+2]); 
%     P = mat_C;
    mat_C = mat2tens(mat_C.',[2^(l-1)+2 2^(l-1)+2 1],3,1:2); 
    
    TTNO{end} = tensor(mat_C,[2^(l-1)+2 2^(l-1)+2 1]);
    
end

end


