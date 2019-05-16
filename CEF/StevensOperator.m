function Bkm = StevensOperator(k,m,Lx,Ly,Lz,L)

cplus   =   1/2; 
cmin    =   1/(2*1i);
Lplus   =   Lx + 1i * Ly;
Lmin    =   Lx - 1i * Ly;
I       =   eye(size(Lx));
X       =   L * (L + 1);

switch k
    case 2
        switch m
            case 2
                Bkm = cplus * (Lplus^2 + Lmin^2);
            case 1
                Bkm = cplus * AntiCommute(Lz, Lplus + Lmin);
            case 0
                Bkm =         (3 * Lz^2 - X * I);
            case -1
                Bkm =  cmin * AntiCommute(Lz, Lplus - Lmin);
            case -2
                Bkm =  cmin * (Lplus^2 - Lmin^2);
            otherwise
                print('check the m value')
        end
    case 4
        switch m
            case 4
                Bkm = cplus * (Lplus^4 + Lmin^4);
            case 3
                Bkm = cplus * AntiCommute(Lz, Lplus^3 + Lmin^3);
            case 2
                Bkm = cplus * AntiCommute(7 * Lz^2 - (X + 5) * I, Lplus^2 + Lmin^2);
            case 1
                Bkm = cplus * AntiCommute(7 * Lz^3 - (3 * X + 1) * Lz, Lplus + Lmin);
            case 0
                Bkm =         (35 * Lz^4 - (30 * X - 25) * Lz^2 + (3 * X^2 - 6 * X) * I);
            case -1
                Bkm =  cmin * AntiCommute(7 * Lz^3 - (3 * X + 1) * Lz, Lplus - Lmin);
            case -2
                Bkm =  cmin * AntiCommute(7 * Lz^2 - (X + 5) * I, Lplus^2 - Lmin^2);
            case -3
                Bkm =  cmin * AntiCommute(Lz, Lplus^3 - Lmin^3);
            case -4
                Bkm =  cmin * (Lplus^4 - Lmin^4);
            otherwise
                print('check the m value')
        end
    otherwise
        print('check the k value')
end
    function C = AntiCommute(A,B)
        C = (1/2) * (A * B + B * A);
    end
end