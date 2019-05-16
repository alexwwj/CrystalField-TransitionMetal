function A1vec = SymmetryAdapetedFunCal(filename,X,Z,k)

dat    = importdata(filename);
cn     = dat.data(:,1);
p(:,1) = dat.data(:,2);
p(:,2) = dat.data(:,3);
p(:,3) = dat.data(:,4);
pm     = dat.data(:,5);
I      = dat.data(:,6);

for ii = 1:length(cn)
    EulerAngle(ii,:) = EulerAngleCalc(Z,X,p(ii,:),cn(ii),pm(ii),1e-12);
    for n = -k:k
        for m = -k:k
            alpha = EulerAngle(ii,1);
            beta  = EulerAngle(ii,2);
            gamma = EulerAngle(ii,3);
            Wigner(n+k+1,m+k+1,ii) = ((-1)^I(ii))^k * WignerD(k,n,m,alpha,beta,gamma);
        end
    end
end

Wmat = sum(Wigner,3)/length(cn);

indx = find(vecnorm(Wmat)>1e-12);

if isempty(indx)
    disp('symmetry adapted functions does not exist')
    A1vec = zeros(2*k+1,1);
else
    A1vec = Wmat(:,indx)./vecnorm(Wmat(:,indx));
end

function WignerD_nm = WignerD(k,n,m,alpha,beta,gamma)
Cnm = (1i)^(abs(n)+n)*(1i)^(-abs(m)-m);

imin = max([0,m-n]);
imax = min([k-n,k+m]);

dk_beta_mn = 0;
for i = imin:imax
    f1 = factorial(k+n) * factorial(k+m) * factorial(k-n) * factorial(k-m);
    f2 = factorial(k-n-i) * factorial(k+m-i) * factorial(i) * factorial(i-m+n); 
    dk_beta_mn = dk_beta_mn + (-1)^(i-m+n) * sqrt(f1)/f2 * ...
        (cos(beta/2))^(2*k+m-n-2*i) * (sin(beta/2))^(2*i+n-m); 
end  

WignerD_nm = Cnm * exp(-1i*n*gamma) * dk_beta_mn * exp(-1i*m*alpha);
end
end