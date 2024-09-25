function [TTNO] = TTNO_no_structure_arbitrary_tree(A,V,X,l,num_l,n,pos)
% This code works for long-range operators of the form
% \sum_i<j A^(i)A^(j)

% note ordering: [Id,a^1,...a1^(m+1),interactions] -> different to binary
% case

I = eye(n(1),n(1));
I = I(:);

m = length(X) - 2;
TTNO = cell(1,m+2);

count = 1;
for kk=1:m
    if iscell(X{kk}) == 1
        count_old = count;
        count = count + count_leaves(X{kk}) - 1;
        TTNO{kk} = TTNO_no_structure_arbitrary_tree(A(count_old:count),V,X{kk},l-1,num_l,n(count_old:count),pos(count_old:count));
        
        count = count + 1;
    else % leaf
        tmp = A{count};
        U = [I tmp(:)];
        TTNO{kk} = U;
        count = count + 1;
    end
end

% connecting tensor
m = length(X) - 2;
s = zeros(1,m+1);
leaves = zeros(m,1);
for kk=1:m
    if iscell(X{kk}) == 1
        s(kk) = count_leaves(X{kk}) + 2;
        leaves(kk) = count_leaves(X{kk});
    else
        s(kk) = 2;
        leaves(kk) = 1;
    end
end

if l == num_l
    % root tensor
    s(m+1) = 1;
    mat_C = zeros(prod(s(1:(end-1))),1);
else
    % other tensors
    s(m+1) = count_leaves(X)+2;
    mat_C = zeros(prod(s(1:(end-1))),s(m+1));
    mat_C(1,1) = 1; % identity
    
    for kk=1:count_leaves(X)
        [tree,num_in_tree] = det_leaf(X,kk);
        prod1 = prod(s(1:(tree-1)));
        
        if tree == 1
            mat_C(num_in_tree,kk+1) = 1;         % A^(i)
        else        
            mat_C(prod1*(num_in_tree-1) + 1,kk+1) = 1; % A^(i)
        end
    end
end

for kk=1:m
    if 1 == iscell(X{kk})
        mat_C(prod(s(1:kk))-prod(s(1:(kk-1))) + 1,end) = 1; % V(i,j)A^(i)A^(j) from the same tree
    end
end


for ii=1:count_leaves(X)
    for jj=1:count_leaves(X)
        [tree1,num_in_tree1] = det_leaf(X,ii);
        [tree2,num_in_tree2] = det_leaf(X,jj);
        if (ii<jj) && (tree1~=tree2)
            % ii-th position "left" tree
            k1 = prod(s(1:(tree1-1)))*(num_in_tree1 - 1) + 1;
            % jj-th position from the "right" tree
            k2 = prod(s((tree1+1):(tree2-1)))*(num_in_tree2 - 1) + 1;
            
            tmp = prod(s(1:tree1))*(k2 - 1) + k1;
            if l == num_l
                mat_C(tmp) = V(pos(ii),pos(jj)); % V(i,j)A^(i)A^(j) from different trees
            else
                mat_C(tmp,end) = V(pos(ii),pos(jj)); % V(i,j)A^(i)A^(j) from different trees
            end
        end
    end
end

TTNO{end-1} = 1;
TTNO{end} = mat2tens(mat_C.',s,m+1,1:m);
TTNO{end} = tensor(TTNO{end},s);

end



function [tree,num_in_tree] = det_leaf(X,n)

m = length(X) - 2;

count = 0;
for ii=1:m
    count_old = count;
    if iscell(X{ii}) == 1
        count = count + count_leaves(X{ii});
    else
        count = count + 1;
    end
    
    if n <= count
        tree = ii;
        num_in_tree = n - count_old + 1;
        break
    end
    
end


end

function [erg] = count_leaves(X)

m = length(X) - 2;
arr = zeros(m+1,1);

for ii=1:m
    if 1==iscell(X{ii})
        tmp = count_leaves(X{ii});
        arr(ii) = tmp(end);
    else
        arr(ii) = 1;
    end
end
arr(end) = sum(arr(1:end-1));
erg = arr(end);

end
