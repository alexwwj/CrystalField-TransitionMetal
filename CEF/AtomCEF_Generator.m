function AtomCEF = AtomCEF_Generator(atomname,L,S,k,m,val,lambda,B,operator)


Basis   =   AngularMomentMatrix(L,S);
Mx      =   Basis.Mx;
My      =   Basis.My;
Mz      =   Basis.Mz;
Jx      =   Basis.Jx;
Jy      =   Basis.Jy;
Jz      =   Basis.Jz;
Jsquare =   Jx^2 + Jy^2 + Jz^2;
N       =   (2 * Basis.LS(1) + 1) * (2 * Basis.LS(2) + 1);
M       =   cat(3,Mx,My,Mz);
J       =   cat(3,Jx,Jy,Jz);

HCEF    =   CrystalFieldHamiltonian(k,m,val,Basis);
HSOC    =   SOCHamiltonian(lambda,Basis);
% heff    =   eps * 5.7883818012e-5 * 1e3;
% Hfield  =   heff * Mx + heff * My + heff * Mz;
muB     =   5.7883818012e-2; % meV / T
Hfield  =   muB * (B(1) * Mx + B(2) * Mx + B(3) * Mx);
H       =   HCEF + HSOC + Hfield; 
tol     =   1e-10; % default value

[CEFval,CEFvec]   =   HamiltonianSolver(H,tol);
% ---- Diagonlize degenerate basis with respect to Operator ----
EvalUnqiue = uniquetol(CEFval,tol);

for i = 1:length(EvalUnqiue)
    indx    =  find(abs(CEFval-EvalUnqiue(i))<tol);
    num     =  length(indx);
    Vec     =  CEFvec(:,indx);
    switch operator
        case 'Mx'
            Op      =  Vec' * (Mx) * Vec;
        case 'My'
            Op      =  Vec' * (My) * Vec;
        case 'Mz'
            Op      =  Vec' * (Mz) * Vec;
        case 'Jx'
            Op      =  Vec' * (Jx) * Vec;
        case 'Jy'
            Op      =  Vec' * (Jy) * Vec;
        case 'Jz'
            Op      =  Vec' * (Jz) * Vec;
        case 'Jsquare'
            Op      =  Vec' * (Jsquare) * Vec;
        case 'I'
            Op      =  Vec' * Vec;
        otherwise 
            disp('Choose an operator from [Mx My Mz] ot [Jx Jy Jz]')
    end
    [~,Mvec]= HamiltonianSolver(Op,tol);
    for j = 1:num
        VecBN(:,indx(j)) = permute(sum(Mvec(:,j) .* permute(Vec,[2,1]),1),[2,1]);
    end
end

AtomCEF.name    = atomname;
AtomCEF.H       = H;
AtomCEF.M       = M;
AtomCEF.J       = J;
AtomCEF.LSdim   = N;
AtomCEF.CEFvec  = gpuArray(VecBN);
AtomCEF.CEFval  = CEFval;
AtomCEF.diffE   = CEFval - min(CEFval);

