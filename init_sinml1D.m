%{
Copyright © 2019 Alexey A. Shcherbakov. All rights reserved.

This file is part of gsmcc.

gsmcc is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

gsmcc is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with sphereml. If not, see <https://www.gnu.org/licenses/>.
%}
%%
function init_sinml1D(kg, ab, ni, ka, kt, esl, es, ec)
    global no;
    global ns;
    global hh;
    global cwz;
    global el;
    global MG;
    
    N = 2^nextpow2(2*no-1);
    
    if ni > 1
        kt = [ka(1,1),kt,ka(1,ni)];% ni+1
        esl = [es, esl, ec];
    else
        kt = [ka(1,1),ka(1,ni)];
        esl = [es, ec];
    end
    ka = [0,ka,0];% ni+1
    
    hh = 0.0;
    for k=1:(ni+1)
        hh = hh + kt(1,k);
    end
    
    tv = 0.5*(hh/ab - hh);
    kt(1,1) = kt(1,1) + tv; kt(1,ni+1) = kt(1,ni+1) + tv;
    
    Z = zeros(1,ni+2);
    Z(1,1) = -0.5*hh;
    for k=2:(ni+2)
        Z(1,k) = Z(1,k-1) + kt(1,k-1);
    end
    
    cwz = zeros(2,ns);
    el = zeros(1,ns);
    
    kk = 1;
    for k=1:ns
        cwz(1,k) = -0.5*hh + hh*(k-0.5)/ns;
        cwz(2,k) = hh/ns;
        if cwz(1,k) > Z(1,kk+1)
            while ((cwz(1,k) > Z(1,kk+1)) && (kk < ni+2))
                kk = kk + 1;
            end
        end
        el(1,k) = esl(1,kk);
    end
    
    MG = zeros(3,ns,N);

    kk = 1;
    for k=1:ns
        z = cwz(1,k); 
        if z > Z(1,kk+1)
            while z > Z(1,kk+1)
                kk = kk + 1;
            end
        end
    
		tt = -0.5*(ka(1,kk+1) - ka(1,kk))/(Z(1,kk+1) - Z(1,kk));
		tv = kg*(ka(1,kk)*(Z(1,kk+1) - z) + ka(1,kk+1)*(z - Z(1,kk)))/(Z(1,kk+1) - Z(1,kk));
		tk = 0.5*tv*tv/(1.0 + 0.5*tv*tv); tn = 2.0/tv;
		MG(1,k,1) = 1.0; MG(1,k,2) = -tt; MG(1,k,N) = -tt;
		t2 = sqrt((1.0-tk)/(1.0+tk)); tm = (1.0-sqrt(1.0-tk*tk))/tk;
		MG(2,k,1) = t2;
        
        tv = tm; m = 1;
        while (2*m < no)
            MG(2,k,2*m+1) = t2*tv; MG(2,k,N-2*m+1) = MG(2,k,2*m+1);
            MG(3,k,2*m+1) = 1i*tn*tt*tv; MG(3,k,N-2*m+1) = -MG(3,k,2*m+1);
			tv = tv*tm; m = m + 1;
        end

		tv = 1.0; m = 0;
        while ((2*m+1) < no)
            MG(2,k,2*m+2) = tt*(tm-1.)*tv; MG(2,k,N-2*m) = MG(2,k,2*m+2);
            MG(3,k,2*m+2) = 0.5*1i*tn*(t2-1.)*tv; MG(3,k,N-2*m) = -MG(3,k,2*m+2);
			tv = tv*tm; m = m + 1;
        end
        
        MG(1,k,:) = fft(MG(1,k,:),N);
        MG(2,k,:) = fft(MG(2,k,:),N);
        MG(3,k,:) = fft(MG(3,k,:),N);
    end
end
