function BolzmanFac = Bolzman(e,temp)

e0 = e-min(e);

nE = length(e0);
nT = length(temp);

E = repmat(e0,1,nT);
T = repmat(temp,nE,1);

kb         =  1/11.604505;
Z          =  repmat(sum(exp(-E./(kb.*T)),1),nE,1);
BolzmanFac =  exp(-E./(kb.*T))./Z;
end