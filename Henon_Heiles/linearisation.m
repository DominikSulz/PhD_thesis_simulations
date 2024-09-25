function [B] = linearisation(D,M1,M2,M3,W,d,lambda)
% This function saves the matrices, which are needed to represent the
% discrete Hamiltonian, in an array B

B = cell(3,d);

for ii=1:d
%     [n,~] = size(D{ii});
%     Z = zeros(n,n);
    % mode k
    if ii==1
        A = D{ii} + 0.5*M1{ii} + W{ii};
    else
        A = D{ii} + 0.5*M1{ii} + W{ii} - (1/3)*lambda*M2{ii-1};
    end
    B{ii,ii} = A;
end

for ii=1:d-1
    
    B{d+ii,ii} = lambda*M1{ii};
    B{d+ii,ii+1} = M3{ii};
end

% old way when using make_operator_kressner
% B = cell(3,d);
% 
% for ii=1:d
%     % first row of H 
%     [n,~] = size(D{ii});
%     Z = zeros(n,n);
%     % mode k
%     if ii==1
%         A = D{ii} + 0.5*M1{ii} + W{ii};
%     else
%         A = D{ii} + 0.5*M1{ii} + W{ii} - (1/3)*lamda*M2{ii-1};
%     end
%     B{1,ii} = A;
%     
%     % second row of H 
%     if ii==d
%         B{2,ii} = Z;
%     else
%         B{2,ii} = lamda*M1{ii};
%     end
%     % third row
%     if ii==1
%         B{3,ii} = Z;
%     else
%         B{3,ii} = M3{ii-1};
%     end
% end
end
