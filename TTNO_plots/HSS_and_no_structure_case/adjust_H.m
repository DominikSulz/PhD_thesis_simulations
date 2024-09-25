function [H,mi_U,mi_V] = adjust_H(H)

% mi = -1 -- U is real but negativ
% mi = -2 -- V is real but negative 
% mi = -3 -- U is complex and negative sign
% mi = -4 -- U is complex and positive sign
% mi = -5 -- V is complex and negative sign
% mi = -6 -- V is complex and positive sign

mi_U = 1;
mi_V = 1;

if H.leafnode == 1
    if H.U < 0 
        H.U = -H.U;
        mi_U = -1;
    end
    if imag(H.U) ~= 0
        if imag(H.U) < 0 
            H.U = 1i*H.U;
            mi_U = -3;
        end
        if imag(H.U) > 0 
            H.U = -1i*H.U;
            mi_U = -4;
        end
    end
    
    if H.V < 0 
        H.V = -H.V;
        mi_V = -2;
    end
    if imag(H.V) ~= 0
        if imag(H.V) < 0 
            H.V = 1i*H.V;
            mi_V = -5;
        end
        if imag(H.U) > 0 
            H.V = -1i*H.V;
            mi_V = -6;
        end
    end
else 
    [H.A11,mi_U,mi_V] = adjust_H(H.A11);
    if mi_U==-1
        H.B12 = -H.B12;
    end
    if mi_U==-3 
        H.B12 = -1i*H.B12;
    end
    if mi_U==-4 
        H.B12 = 1i*H.B12;
    end
    
    if mi_V==-2
        H.B21 = -H.B21;
    end
    if mi_V==-5 
        H.B21 = 1i*H.B21; % old H.B21 = -1i*H.B21;
    end
    if mi_V==-6 
        H.B21 = -1i*H.B21; % old H.B21 = 1i*H.B21;
    end
    
    [H.A22,mi_U,mi_V] = adjust_H(H.A22);
    if mi_U==-1
        H.B12 = -H.B12;
    end
    if mi_U==-3 
        H.B12 = -1i*H.B12;
    end
    if mi_U==-4 
        H.B12 = 1i*H.B12;
    end
    
    if mi_V==-2
        H.B21 = -H.B21;
    end
    if mi_V==-5 
        H.B21 = 1i*H.B21; % old H.B21 = -1i*H.B21;
    end
    if mi_V==-6 
        H.B21 = -1i*H.B21; % old H.B21 = 1i*H.B21;
    end
end

end

%% old
% if H.leafnode == 1
%     mi = 1;
%     if H.U < 0 
%         H.U = -H.U;
%         mi = mi*(-1);
%     end
% else 
%     [H.A11,mi] = adjust_H(H.A11);
%     if mi==-1
%         H.B12 = -H.B12;
%     end
%     [H.A22,mi] = adjust_H(H.A22);
%     if mi==-1
%         H.B21 = -H.B21;
%     end
% end


%% old 2
% if H.leafnode == 1
%     mi = 1;
%     if (H.U < 0) 
%         H.U = -H.U;
%         mi = mi*(-1);
%     end
%     if imag(H.U) < 0
%         if imag(H.V) > 0
%             H.U = -H.U;
%             mi = mi*(-2);
%         end
%     end
%     if imag(H.U) > 0
%         if  imag(H.V) < 0
%             H.V = -H.V;
%             mi = mi*(-2);
%         end
%     end
% else 
%     [H.A11,mi] = adjust_H(H.A11);
%     if mi==-1
%         H.B12 = -H.B12;
%     end
%     if mi==-2
%         H.B12 = -H.B12;
%     end
%     [H.A22,mi] = adjust_H(H.A22);
%     if mi==-1
%         H.B21 = -H.B21;
%     end
%     if mi==-2
%         H.B21 = -H.B21;
%     end
% end