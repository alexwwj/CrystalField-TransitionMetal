clear
clc 

atomname =   'Co2';
L        =   3; 
S        =   3/2;

k        =   [2,4];
m        =   {[0],[0,3]};
val      =   {[25.85],[0.4555,14.017]};
lambda   =   -1.021 * 22.32;

T        =   200;
Ei       =   1e2;
theta    =   [3,134];
Evec     =   linspace(0,Ei,5e2);
Qvec     =   linspace(0,8,5e2);
num      =   2e3;
FWHMp    =   [0 0 -0.01 12];

AtomCEF  =   AtomCEF_Generator(atomname,L,S,k,m,val,lambda);
PowderSpec = PowderSpectrum(AtomCEF,T,Ei,theta,Qvec,Evec,num,FWHMp);
%%
figure(1)
pcolor(PowderSpec.Qvec,PowderSpec.Evec,PowderSpec.spec)
shading flat;
set(gca,'YDir','normal')
% title(sprintf('Concoluted powder spectra: Re Sperp($\frac{1}{2}Omega$,Q), T = %4.2f',T),'interpreter','latex')
colorbar
caxis([0 0.02])
colormap jet
%%
indx = find(0<Qvec<10);
PowderSpec.spec(isnan(PowderSpec.spec)) = 0;
figure(2)
plot(Evec,sum(PowderSpec.spec(:,indx),2))
%%
CEFvec = AtomCEF.CEFvec;
CEFval = AtomCEF.CEFval;
conjv1 = gpuArray(repmat(CEFvec',1,1,3));
v1     = gpuArray(repmat(CEFvec',1,1,3));
M      = gpuArray(AtomCEF.M);


