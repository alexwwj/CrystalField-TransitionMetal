function U = SOCHamiltonian(lambda,Basis)

Lx = Basis.Lx;
Ly = Basis.Ly;
Lz = Basis.Lz;
Sx = Basis.Sx;
Sy = Basis.Sy;
Sz = Basis.Sz;

U = lambda * (Lx * Sx + Ly * Sy + Lz * Sz);
end