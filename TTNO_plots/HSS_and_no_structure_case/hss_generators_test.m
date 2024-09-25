function [U,V] = hss_generators_test(A)
if A.leafnode == 1
    U = A.U;
    V = A.V;
else
    [Ul,Vl] = hss_generators_test(A.A11);
    [Ur,Vr] = hss_generators_test(A.A22);
    U = [ Ul * A.Rl ; Ur * A.Rr ]
    V = [ Vl * A.Wl ; Vr * A.Wr ]
end
end