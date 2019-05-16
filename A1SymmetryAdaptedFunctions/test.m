clear
clc
filename = 'pg_D3d.dat';
X = [1 0 0];
Z = [0 0 1];
k = 4;
A1vec = SymmetryAdapetedFunCal(filename,X,Z,k)
