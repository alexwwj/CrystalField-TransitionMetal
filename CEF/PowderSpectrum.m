function NeutronPowderSpec = PowderSpectrum(AtomCEF,T,Ei,theta,Qvec,Evec,num,p)

if ismember(0,Qvec) == 1
    ind = Qvec == 0;
    Qvec(ind) = eps;
end

nQ = length(Qvec);
nE = length(Evec);
PowderSpec = zeros(nE,nQ);
KinematicSpec = KinematicConstraint(Ei,theta,Evec,Qvec);

% ------ Randam Q generation ------
TH = 2 * pi * rand(nQ,num);
PH = asin(-1 + 2 * rand(nQ,num));
[Qx,Qy,Qz] = sph2cart(TH,PH,Qvec');

% ------ Neutron powder spectrum ------
h = waitbar(0,'Please wait...');
for i = 1:nQ
    waitbar(i / nQ,h,sprintf('Already complete %3.2f%%, please wait...',i/nQ*100))
    Q = [Qx(i,:);Qy(i,:);Qz(i,:)];
    PowderSpec(:,i) = PowderCrossSection(AtomCEF,T,Ei,Q,Evec,p);
end
NeutronPowderSpec.spec  =  PowderSpec .* KinematicSpec;
NeutronPowderSpec.Qvec  =  Qvec;
NeutronPowderSpec.Evec  =  Evec;
NeutronPowderSpec.theta =  theta;
NeutronPowderSpec.T     =  T;
NeutronPowderSpec.Ei    =  Ei;

close(h)

end

function spec = KinematicConstraint(Ei,theta,Evec,Qvec)

spec    = double(bsxfun(@(Q,E)fun(Q,E,Ei,theta),Qvec,Evec'));
spec(spec==0) = NaN;

function c = fun(Q,Evec,Ei,theta)
mu      =   1.660539040e-27; 
m       =   1.008644904 * mu;
h       =   6.62607004e-34; 
const   =   (h / pi / 2) ^2 / (2 * m);
ev      =   1.6021766208e-19;
ei      =   Ei * 1e-3 * ev; 
E       =   Evec * 1e-3 * ev;

f1      =   2 * ei - E - 2 * (ei * (ei - E)) .^ (1/2) * cosd(theta(1)) - const * (Q * 1e10) .^2;
f2      =   2 * ei - E - 2 * (ei * (ei - E)) .^ (1/2) * cosd(theta(end)) - const * (Q * 1e10) .^2;
c       =   (f1 .* f2)<0;
end

end

