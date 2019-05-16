clear
clc 

%%%% ---- Input single ion parameters ---- %%%%
atomname =   'Co2';
L        =   3; 
S        =   3/2;

k        =   [2,4];
m        =   {[0],[0,3]};
val      =   {[25.85],[0.4555,14.017]};
lambda   =   -1.021 * 22.32;
B        =   [0 0 0];

T  =  linspace(1,300,1e3);
kb =  1/11.604505;
%%%% ---- calculated ground state multiplet splittings ---- %%%%
AtomCEF  =   AtomCEF_Generator(atomname,L,S,k,m,val,lambda,B,'Jz');

%%%% ---- DC susceptibility calculation ---- %%%%
CEFvec = gather(AtomCEF.CEFvec);
CEFval = AtomCEF.CEFval;
m      = gpuArray(repmat(CEFvec',1,1,3));
n      = gpuArray(repmat(CEFvec,1,1,3));
M      = gpuArray(AtomCEF.M);
N      = AtomCEF.LSdim;

%%%% ---- <m|L+2S|n><n|L+2S|m> ---- %%%%
Mmat = gather(pagefun(@mtimes,m,pagefun(@mtimes,M,n)));
mx   = Mmat(:,:,1);
my   = Mmat(:,:,2);
mz   = Mmat(:,:,3);

Magxx = mx .* conj(mx);
Magxy = mx .* conj(my);
Magxz = mx .* conj(mz);

Magyx = my .* conj(mx);
Magyy = my .* conj(my);
Magyz = my .* conj(mz);

Magzx = mz .* conj(mx);
Magzy = mz .* conj(my);
Magzz = mz .* conj(mz);

Mag = cat(3,Magxx,Magxy,Magxz,Magyx,Magyy,Magyz,Magzx,Magzy,Magzz);

%%%% ---- Create masks to distinguish the Van-veleck and dynamical ---- %%%%
%%%% ---- contributions                                            ---- %%%%
Em      = repmat(CEFval,1,N);
En      = repmat(CEFval',N,1);
deltaE  = Em-En;
% deltaE(abs(deltaE)<1e-12) = 0;

mask_VanVleck   = double(abs(deltaE)<1e-12);
mask_Dynamical  = ones(N,N)-mask_VanVleck;

nT = length(T);
chi_VanVleck    = repmat(mask_VanVleck .* Mag,1,1,1,nT);
chi_Curie       = repmat(mask_Dynamical .* Mag,1,1,1,nT);

%%%% ---- Calculate the Van-Vleck susceptibility part ---- %%%%
BolzmanFac      =  reshape(Bolzman(CEFval,T),N,1,1,nT);
Pn_VanVleck     =  repmat(BolzmanFac,1,N,9,1);
Chi_VanVleck1   =  reshape(sum(reshape(Pn_VanVleck .* chi_VanVleck,N*N,9,nT),1),9,nT)./repmat(T,9,1)/kb;

MeanMx = sum(repmat(diag(mx),1,nT) .* Bolzman(CEFval,T),1);
MeanMy = sum(repmat(diag(my),1,nT) .* Bolzman(CEFval,T),1);
MeanMz = sum(repmat(diag(mz),1,nT) .* Bolzman(CEFval,T),1);

Chi_VanVleck2   =  [   MeanMx .* MeanMx ./ T; MeanMx .* MeanMy ./ T; MeanMx .* MeanMz ./ T; 
                       MeanMy .* MeanMx ./ T; MeanMy .* MeanMy ./ T; MeanMy .* MeanMz ./ T;
                       MeanMz .* MeanMx ./ T; MeanMz .* MeanMy ./ T; MeanMz .* MeanMz ./ T;   ]/kb;
                 
Chi_VanVleck    =  Chi_VanVleck1 - Chi_VanVleck2;                

%%%% ---- Calculate the Curie susceptibility part ---- %%%%
Pm_Curie    =  repmat(BolzmanFac,1,N,9,1);
Pn_Curie    =  permute(Pm_Curie,[2,1,3,4]);
Pmn_Curie   =  Pm_Curie-Pn_Curie;
dEnm_Curie  =  repmat(-deltaE,1,1,9,nT);
dEnm_Curie(abs(dEnm_Curie)<1e-12)  =  eps;
Chi_Curie   =  reshape(sum(reshape(Pmn_Curie .* chi_Curie ./ dEnm_Curie,N*N,9,nT),1),9,nT);

muB      =  9.274009994e-24; % J * T^-1
NA       =  6.02214076e23;  
mev2joul =  1.60218e-22;
mu0      =  4*pi*1e-7; % T * m/A
factor   =  NA * muB^2 * mu0 / mev2joul / (4*pi*1e-6); %emu mol-1
% factor = 1;
Chi_ion  =  permute(reshape(factor * real(Chi_Curie + Chi_VanVleck),3,3,nT),[2,1,3]);
% Chi_ion(Chi_ion<1e-12) = 0;

ChiPow = reshape((Chi_ion(1,1,:) + Chi_ion(2,2,:) + Chi_ion(3,3,:))/3,1,nT);
lambda = 23.99;
ChiMF = (-3.1e-4 + ChiPow ./ (1+lambda*ChiPow));


figure(1)
plot(T,1./ChiMF)


