function Q = GramSchmidt(A)
[m,n] = size(A);
   R  = zeros(m,n);
   Q  = zeros(m,n);
for j = 1:n
    v = A(:,j);
    for i = 1:j-1
        R(i,j) = Q(:,i)' * A(:,j);
        v      = v - R(i,j) * Q(:,i);
    end
    R(j,j) = norm(v);
    Q(:,j) = v / R(j,j);
end