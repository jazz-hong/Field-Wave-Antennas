%**************************************************************************
%  Diffraction coefficients caculation function
%  Based on the Fortran program in Balanis' book, subrouting wdc().
%  Some modifications have been made due to Matlab issues.
%  If you find any bugs in the Matlab program or you have suggestions,
%       please contact me at zyun at hawaii dot edu.
%  
%  Zhengqing Yun, 2000~2006
%**************************************************************************
function [dcs,dch,D] = wdcpub(R,phid,phipd,btd,fn)
sml = 0.001;
if abs(fn-0.5) < sml; return; end

d2r = pi/180.0;
ct = -exp(-pi/4*j)/(4*pi); 
tpi = 2*pi;	utpi = 1.0/tpi;	

phr = phid*d2r; phpr = phipd*d2r; sbo = sin(d2r*btd);
dkl = tpi*R; ufn = 1./fn; betar = phr - phpr; term = ct/fn/sbo;

D = []; ft = [];
sgn = 1.0;
for i=1:4;
    N = round(ufn*(0.5*sgn + utpi*betar));
    ang = tpi*fn*N - betar;
    g = 1.0+cos(ang);
    X = dkl*g; Y = pi + sgn*betar;

    if abs(X) < 1.0e-10
        ft(i) = 0.0;
        if abs(cos(0.5*Y*ufn)) >= sml;
            sny = sign(sign(Y)-0.5);
            ft(i) = sny*sqrt(dkl)*(1.7725+1.7725*j);
            ft(i) = (ft(i) - (0.0+2.0*j)*dkl*(pi-sgn*ang))*fn;
        end
    else
        ft(i) = ftfpub(X)/tan(0.5*Y*ufn);
    end

    sgn = -sgn;
    D(i) = term*ft(i);
    if (sgn >= 0.0)
        betar = phr + phpr;
    end
end
dcs = D(1)+D(2)-D(3)-D(4); %soft polarization
dch = D(1)+D(2)+D(3)+D(4); %hard polarization
