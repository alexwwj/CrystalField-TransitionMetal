clear
clc
close all
atomname =   'Co2';
L        =   3; 
S        =   3/2;

k        =   [2,4];
m        =   {[0],[0,3]};
val      =   {[25.85],[0.4555,14.017]};
lambda   =   -1.021 * 22.32;
B        =   [0 0 0];

AtomCEF  =   AtomCEF_Generator(atomname,L,S,k,m,val,lambda,B,'Jz');

T        =   linspace(1,300,1e3);
Chi_dia  =   -3.1e-4;
lambda   =   23.99;
Chi      =   SingleIonMagSusceptibility(AtomCEF,T,Chi_dia,lambda);

figure(1)
plot(T,1./Chi.ChiPowMF)
xlabel('Temperature [K]')
ylabel('1/\chi [mol Oe/emu]')
%%
T        =   5;
Ei       =   250;
theta    =   [3,134];
Evec     =   linspace(0,Ei,5e2);
Qvec     =   linspace(0,8,5e2);
num      =   1e3;
FWHMp    =   [0 0 0 13];

PowderSpec = PowderSpectrum(AtomCEF,T,Ei,theta,Qvec,Evec,num,FWHMp);
%%
figure(2)
pcolor(PowderSpec.Qvec,PowderSpec.Evec,PowderSpec.spec)
shading flat;
set(gca,'YDir','normal')
% title(sprintf('Concoluted powder spectra: Re Sperp($\frac{1}{2}Omega$,Q), T = %4.2f',T),'interpreter','latex')
colorbar
caxis([0 0.08])
colormap jet

%%
qmin = 4;
qmax = 5;
indx = find((qmin<Qvec) & (Qvec<qmax));
% PowderSpec.spec(isnan(PowderSpec.spec)) = 0;
figure(3)
plot(Evec,nanmean(PowderSpec.spec(:,indx),2))

% figure(2)
% imagesc(PowderSpec.Qvec,PowderSpec.Evec,PowderSpec.spec)