clear
clc
close all

Z = [0 0 1];
X = [1 0 0];

p = [-1 -1 1];
n = 3;
pm = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = Z/norm(Z);
X = X/norm(X);
Y = cross(Z,X);

Xp = dot(p,X);
Yp = dot(p,Y);
Zp = dot(p,Z);
u = [Xp Yp Zp]/norm([Xp Yp Zp]);

theta  = (-1)^(pm) * 2 * pi / n;

R(1,1) = cos(theta) + u(1)^2 * (1-cos(theta));
R(1,2) = u(1) * u(2) * (1-cos(theta)) - u(3) * sin(theta);
R(1,3) = u(1) * u(3) * (1-cos(theta)) + u(2) * sin(theta);

R(2,1) = u(2) * u(1) * (1-cos(theta)) + u(3) * sin(theta);
R(2,2) =  cos(theta) + u(2)^2 * (1-cos(theta));
R(2,3) = u(2) * u(3) * (1-cos(theta)) - u(1) * sin(theta);

R(3,1) = u(3) * u(1) * (1-cos(theta)) - u(2) * sin(theta); 
R(3,2) = u(3) * u(2) * (1-cos(theta)) + u(1) * sin(theta);
R(3,3) = cos(theta) + u(3)^2 * (1-cos(theta));


if R(3,3)<1
    if R(3,3)>-1
        thetaY = acos(R(3,3));
        thetaZ0 = atan2(R(2,3),R(1,3));
        thetaZ1 = atan2(R(3,2),-R(3,1));
    else
        %%% not a unique solution thetaZ1 - thetaZ0 = atan2(R21,R22)
        thetaY = pi;
        thetaZ1 = atan2(R(2,1),R(2,2));
        thetaZ0 = 0;
    end
else
    %%% not a unique solution thetaZ1 + thetaZ0 = atan2(R21,R22) 
    thetaY = 0;
    thetaZ0 = 0;
    thetaZ1 = atan2(R(2,1),R(2,2));
end

EulerAngle(1) = thetaZ1 * (thetaZ1>=0 && thetaZ1<2*pi) + (2 * pi + thetaZ1) * (thetaZ1<0) + (-2 * pi + thetaZ1) * (thetaZ1==2*pi);
EulerAngle(2) = thetaY;
EulerAngle(3) = thetaZ0 * (thetaZ0>=0 && thetaZ0<2*pi) + (2 * pi + thetaZ0) * (thetaZ0<0) + (-2 * pi + thetaZ0) * (thetaZ0==2*pi);
