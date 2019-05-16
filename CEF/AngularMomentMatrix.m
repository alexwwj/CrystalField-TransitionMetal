function Basis = AngularMomentMatrix(L,S)

nL = 2 * L + 1;
nS = 2 * S + 1;
N  = nL * nS;

for i = 1:nL
    for ii = 1:nS
        phi{i,ii} = [-L+i-1,-S+ii-1];
    end
end

phi = reshape(phi,[N,1]);

Lx = zeros(N,N);
Ly = zeros(N,N);
Lz = zeros(N,N);
Sx = zeros(N,N);
Sy = zeros(N,N);
Sz = zeros(N,N);

for i = 1:N
    Lzi = phi{i}(1);
    Szi = phi{i}(2);
    for j = 1:N
        Lzj = phi{j}(1);
        Szj = phi{j}(2);
        Lx(i,j) = ((Delta(Lzi,Lzj + 1) + Delta(Lzi + 1,Lzj)) * (1 / 2) * sqrt(L * (L + 1) - Lzi * Lzj)) * Delta(Szi,Szj);
        Ly(i,j) = ((Delta(Lzi,Lzj + 1) - Delta(Lzi + 1,Lzj)) * (1 / (2 * 1i)) * sqrt(L * (L + 1) - Lzi * Lzj)) * Delta(Szi,Szj);
        Lz(i,j) = (Delta(Lzi,Lzj) * Lzi) * Delta(Szi,Szj);
        Sx(i,j) = ((Delta(Szi,Szj + 1) + Delta(Szi + 1,Szj)) * (1 / 2) * sqrt(S * (S + 1) - Szi * Szj)) * Delta(Lzi,Lzj);
        Sy(i,j) = ((Delta(Szi,Szj + 1) - Delta(Szi + 1,Szj)) * (1 / (2 * 1i)) * sqrt(S * (S + 1) - Szi * Szj)) * Delta(Lzi,Lzj);
        Sz(i,j) = (Delta(Szi,Szj) * Szi) * Delta(Lzi,Lzj);
    end
end

Basis.phi =  phi;
Basis.LS  =  [L,S];
Basis.Lx  =  Lx;
Basis.Ly  =  Ly;
Basis.Lz  =  Lz;
Basis.Sx  =  Sx;
Basis.Sy  =  Sy;
Basis.Sz  =  Sz;

Mx = Lx + 2* Sx;
My = Ly + 2* Sy;
Mz = Lz + 2* Sz;

Jx = Lx + Sx;
Jy = Ly + Sy;
Jz = Lz + Sz;

Basis.Mx  =  Mx;
Basis.My  =  My;
Basis.Mz  =  Mz;

Basis.Jx  =  Jx;
Basis.Jy  =  Jy;
Basis.Jz  =  Jz;

function value = Delta(m,n)
    if m == n
        value = 1;
    else
        value = 0;
    end
end
end