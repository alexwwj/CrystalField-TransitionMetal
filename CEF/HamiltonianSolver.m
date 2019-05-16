function [E,Q] = HamiltonianSolver(H,tol)

[V,~]   =   eig(H);
Q = GramSchmidt(V);

if (norm(Q'*Q-eye(size(H)))) > tol
    disp('The Gram-Schmidt Othornormalization does not work or the tolerence is too samll')
    return;
end

E = real(diag(Q'*H*Q));

if norm((Q'*H*Q)-diag(E)) > 10*tol
    disp('The othornormal matrix cannot diagonalize the Hamiltonian')
    return;
end

[E,indx] = sort(E);
 Q = Q(:,indx');

end