function [Y1] = F_Hamiltonian(t,Y,A,d)

Y1 = apply_operator_nonglobal(Y,A,d);

end