clear
clc
close all

E_Eulerangle    =     [
                        0 * pi    , 0 * pi    , 0 * pi    ;
                      ];

C3_Eulerangle   =     [   
                        1 * pi / 2, 1 * pi / 2, 0 * pi    ;
                        3 * pi / 2, 1 * pi / 2, 1 * pi    ;
                        1 * pi / 2, 1 * pi / 2, 1 * pi    ;
                        3 * pi / 2, 1 * pi / 2, 0 * pi    ;
                        1 * pi    , 1 * pi / 2, 1 * pi / 2;
                        0 * pi    , 1 * pi / 2, 3 * pi / 2;
                        0 * pi    , 1 * pi / 2, 1 * pi / 2;
                        1 * pi    , 1 * pi / 2, 3 * pi / 2;
                      ];



C4_Eulerangle   =     [   
                        1 * pi / 2, 1 * pi / 2, 3 * pi / 2;
                        0 * pi / 2, 1 * pi / 2, 0 * pi    ;
                        1 * pi / 2, 0 * pi    , 0 * pi    ;
                        3 * pi / 2, 1 * pi / 2, 1 * pi / 2;
                        1 * pi    , 1 * pi / 2, 1 * pi    ;
                        3 * pi / 2, 0 * pi    , 0 * pi    ;
                      ];           
                  
C2abc_Eulerangle   =   [   
                        1 * pi / 2, 1 * pi    , 0 * pi    ;
                        3 * pi / 2, 1 * pi    , 0 * pi    ;
                        1 * pi    , 1 * pi / 2, 0 * pi    ;
                        1 * pi / 2, 1 * pi / 2, 1 * pi / 2;
                        0 * pi    , 1 * pi / 2, 1 * pi    ;
                        3 * pi / 2, 1 * pi / 2, 3 * pi / 2;
                      ];      
                  
                  
C2xyz_Eulerangle   =   [   
                        1 * pi    , 1 * pi    , 0 * pi    ;
                        0 * pi    , 1 * pi    , 0 * pi    ;
                        1 * pi    , 0 * pi    , 0 * pi    
                      ];                        

k = 4;
h = 24;
%%
clear D

WE = zeros(2*k+1,2*k+1);

for i = 1:1
    alpha = E_Eulerangle(i,1);
    beta  = E_Eulerangle(i,2);
    gamma = E_Eulerangle(i,3);
    for n = -k:k
        for m = -k:k
            D{i}(n+k+1,m+k+1) = (1/h) * WignerD(k,n,m,alpha,beta,gamma);
        end
    end
    WE = WE+D{i}; 
end
%%
clear D

WC3 = zeros(2*k+1,2*k+1);

for i = 1:8
    alpha = C3_Eulerangle(i,1);
    beta  = C3_Eulerangle(i,2);
    gamma = C3_Eulerangle(i,3);
    for n = -k:k
        for m = -k:k
            D{i}(n+k+1,m+k+1) = (1/h) * WignerD(k,n,m,alpha,beta,gamma);
        end
    end
    WC3 = WC3+D{i}; 
end
%%
clear D

WC4 = zeros(2*k+1,2*k+1);

for i = 1:6
    alpha = C4_Eulerangle(i,1);
    beta  = C4_Eulerangle(i,2);
    gamma = C4_Eulerangle(i,3);
    for n = -k:k
        for m = -k:k
            D{i}(n+k+1,m+k+1) = (1/h) * WignerD(k,n,m,alpha,beta,gamma);
        end
    end
    WC4 = WC4+D{i}; 
end

%%
clear D

WC2abc = zeros(2*k+1,2*k+1);

for i = 1:6
    alpha = C2abc_Eulerangle(i,1);
    beta  = C2abc_Eulerangle(i,2);
    gamma = C2abc_Eulerangle(i,3);
    for n = -k:k
        for m = -k:k
            D{i}(n+k+1,m+k+1) = (1/h) * WignerD(k,n,m,alpha,beta,gamma);
        end
    end
    WC2abc = WC2abc+D{i}; 
end

%%
clear D

WC2xyz = zeros(2*k+1,2*k+1);

for i = 1:3
    alpha = C2xyz_Eulerangle(i,1);
    beta  = C2xyz_Eulerangle(i,2);
    gamma = C2xyz_Eulerangle(i,3);
    for n = -k:k
        for m = -k:k
            D{i}(n+k+1,m+k+1) = (1/h) * WignerD(k,n,m,alpha,beta,gamma);
        end
    end
    WC2xyz = WC2xyz+D{i}; 
end


% alpha = pi/2;
% beta  = pi/2;
% gamma = pi/2;
% 
% for n = -k:k
%     for m = -k:k
%         Dk(n+k+1,m+k+1) = WignerD(k,n,m,alpha,beta,gamma);
%     end
% end
% 
% print


W = WE+WC3+WC4+WC2abc+WC2xyz
%(W./vecnorm(W))
real(W(:,1)/vecnorm(W(:,1)))



















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