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
% Ei       =   1e2;
% theta    =   [3,134];
% Evec     =   linspace(0,Ei,5e3);
% Qvec     =   linspace(0,8,5e2);
% num      =   1e3;
% FWHMp    =   [0 0 -0.01 12];

AtomCEF  =   AtomCEF_Generator(atomname,L,S,k,m,val,lambda);

CEFvec = AtomCEF.CEFvec;
CEFval = AtomCEF.CEFval;
pi     = gpuArray(repmat(CEFvec',1,1,3));
pj     = gpuArray(repmat(CEFvec,1,1,3));
M      = gpuArray(AtomCEF.M);
N      = AtomCEF.LSdim;

%%%% ---- f = <pi|L+2S|pj><pj|L+2S|pi> ---- %%%%
f1 = gather(pagefun(@mtimes,pi,pagefun(@mtimes,M,pj)));
f = f1 .* conj(f1); 
Epi = repmat(CEFval,1,N);
Epj = repmat(CEFval',N,1);
deltaE = Epi-Epj;
deltaE(abs(deltaE)<1e-12) = 0;

mask1 = double(deltaE==0);
mask2 = ones(N,N)-mask1;

%%%% ---- chi1 = <pi|L+2S|pj><pj|L+2S|pi> where Epi = Epj---- %%%%
%%%% ---- chi2 = <pi|L+2S|pj><pj|L+2S|pi> where Epi ~= Epj---- %%%%
chi1 = mask1 .* f;
chi2 = mask2 .* f;

T  =  linspace(1e-1,300,1e3);
nT = length(T);
kb =  11.604505;
BolzmanFac = reshape(Bolzman(CEFval,T),N,1,nT);
Pn1 = repmat(BolzmanFac,1,N,1);

Chixx1 = sum(reshape(Pn1 .* chi1(:,:,1),N*N,nT),1)./T/kb;
Chiyy1 = sum(reshape(Pn1 .* chi1(:,:,2),N*N,nT),1)./T/kb;
Chizz1 = sum(reshape(Pn1 .* chi1(:,:,3),N*N,nT),1)./T/kb;

Pnj = repmat(BolzmanFac,1,N,1);
Pni = permute(Pnj,[2,1,3]);
Pnji = Pnj-Pni;

dE = -deltaE;
dE(dE==0)=1e3;

Chixx2 = sum(reshape(Pnji .* chi2(:,:,1) ./ dE,N*N,nT),1);
Chiyy2 = sum(reshape(Pnji .* chi2(:,:,2) ./ dE,N*N,nT),1);
Chizz2 = sum(reshape(Pnji .* chi2(:,:,3) ./ dE,N*N,nT),1);

% figure(1)
% plot(T,Chixx1+Chiyy1+Chizz1)
% 
% figure(2)
% plot(T,Chixx2+Chiyy2+Chizz2)

Chixy_ion = (Chixx1+Chixx2);
lambda = 6;
ChixyMF = Chixy_ion./(1+lambda*Chixy_ion);

Chiz_ion = (Chizz1+Chizz2);
ChizMF = Chiz_ion./(1+lambda*Chiz_ion);


figure(3)
clf; hold on
plot(T,1./ChixyMF)
plot(T,1./ChizMF)
hold off

function BolzmanFac = Bolzman(e,temp)

nE = length(e);
nT = length(temp);

E = repmat(e,1,nT);
T = repmat(temp,nE,1);

kb         =  11.604505;
Z          =  repmat(sum(exp(-E./kb./T),1),nE,1);
BolzmanFac =  exp(-E./kb./T)./Z;
end
