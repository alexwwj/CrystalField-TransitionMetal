function SpecPowd = SqwPowder(atomname,T,Ei,Q,Omega,CEFvec,CEFval,Basis,p)

%%%% ------------- input structure -------------
%%%% atomname   A string, contains the name of the magnetic ion in FillProf
%%%%            notaion.
%%%% T          Temperature at which the scattering is conducted.
%%%% Ei         The incident neutron energy in meV unit.
%%%% Q          Momentum transfers in Angstrom^-1 units with dimension of
%%%%            [3 nQ].
%%%% omega      Energy transfer in meV units with dimension of [1 nomega].
%%%% CEFvec     The solved eigenvectors for the ground state multiplet
%%%%            splittings.
%%%% CEFval     The solved eigenenergies for the ground state multiplets.
%%%% Basis      A structure, contains [Mx,My,Mz] matrices.
%%%% gamma      The width for the Lorentz resolution function. 

% ------ The incident and the scattered momentum transfers ------
[~,nQ]  =  size(Q);
ki      =  [0; 0; sqrt(Ei/2.072)];
kf      =  repmat(ki,1,nQ) + Q;

% ------ Magnetic form factor ------
FormFacVal = sw_mff(atomname,Q);

% ------ The perpendicular magnetic moment components ------
Mx = Basis.Mx;
My = Basis.My;
Mz = Basis.Mz;

N   = (2 * Basis.LS(1) + 1) * (2 * Basis.LS(2) + 1);
M   = cat(3,Mx,My,Mz);
M0  = reshape(M,N*N,3);
M1  = (M0 * Q) ./ vecnorm(Q).^2;
M2  = reshape(M1,N,N,nQ);
Mperp = zeros(N,N,3,nQ);
for i = 1:nQ
    M3 = cat(3,M2(:,:,i)*Q(1,i),M2(:,:,i)*Q(2,i),M2(:,:,i)*Q(3,i));
    Mperp(:,:,:,i) = M - M3;
end

% ------ The energy transfer dependent part of the dynamical structure factor------
nE          = length(Omega);
Evec1       = repmat(CEFval,1,N);
Evec2       = repmat(CEFval',N,1);
DeltaE      = repmat(Evec1-Evec2,1,1,nE);
Evec        = repmat(reshape(Omega,1,1,nE),N,N,1);
EnergyReso  = ResoFun(Evec,DeltaE,p);

% ------ The momentum transfer dependent part of the dynamical structure factor------
gamma       = 1.91;
r           = 2.818e-15; % unit: m
fac         = (gamma * r / 1e-14)^2 * (vecnorm(kf) / vecnorm(ki)) .* abs(FormFacVal / 2).^2;
BolzmanFac  = repmat(Bolzman(CEFval,CEFval,T),1,N,3,nQ); % dim [N N 3 nQ]
mfac        = repmat(reshape(fac,1,1,1,nQ),N,N,3,1);     % dim [N N 3 nQ]
EigVec1     = gpuArray(repmat(CEFvec',1,1,3,nQ));        % dim [N N 3 nQ]
EigVec2     = gpuArray(repmat(CEFvec,1,1,3,nQ));         % dim [N N 3 nQ]
A           = gpuArray(Mperp);
% conjA       = gpuArray(conj(permute(Mperp,[2,1,3,4])));
W           = gather(pagefun(@mtimes,EigVec1,pagefun(@mtimes,A,EigVec2)));
% conjW       = gather(pagefun(@mtimes,EigVec1,pagefun(@mtimes,conjA,EigVec2)));
mfun        = BolzmanFac .* mfac .* W .* conj(W);%permute(conjW,[2,1,3,4]);
sumQ        = reshape(sum(reshape(mfun,N*N,3*nQ),2)/nQ,N,N);

% ------ The double differential neutron cross section for a powder sample------
sumQE       = repmat(sumQ,1,1,nE);
SpecPowd    = real(sum(reshape(sumQE.*EnergyReso,N*N,nE),1));

function fun = ResoFun(x,omega,p)
Gamma =  p(1) * omega.^3 + p(2) * omega.^2 + p(3) * omega.^1 + p(4); 
fun   =  (1/pi) * Gamma ./ (Gamma.^2 + (omega - x).^2);
end

function BolzmanFac = Bolzman(E0,E,T)
kb = 11.604505;
Z = sum(exp(-E/kb/T),1);
BolzmanFac = exp(-E0/kb/T)/Z;
end

end
