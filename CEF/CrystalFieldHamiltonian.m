function U = CrystalFieldHamiltonian(k,m,val,Basis)

%%%% k   = [2,4]
%%%% m   = {[0,1,2],[2,3,4]}
%%%% val = {[10,20,5],[20,30,2]}
Lx = Basis.Lx;
Ly = Basis.Ly;
Lz = Basis.Lz;
L  = Basis.LS(1);

len  =  length(k);
k    =  reshape(k,[1,len]);
m    =  reshape(m,[1,len]);
val  =  reshape(val,[1,len]);

if length(m) ~= len || length(val) ~= len
    disp('check the input CF order parameters')
    return;
end

for i = 1:len
    ki     = k(i);
    lenmi  = length(m{i});
    Bkm{i} = zeros(size(Lz)); 
    for ii = 1:lenmi
        mi   = m{i}(ii);
        vali = val{i}(ii);
        Bkm{i} = Bkm{i} + vali * StevensOperator(ki,mi,Lx,Ly,Lz,L);
    end
end
U = sum(cat(3,Bkm{:}),3);