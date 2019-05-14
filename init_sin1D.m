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
function init_sin1D(kg, ab, eps1, eps2)
    global no;
    global ns;
    global hh;
    global cwz;
    global el;
    global MG;
    
    N = 2^nextpow2(2*no-1);
    
    cwz = zeros(2,ns);
    el = zeros(1,ns);
    for k=1:ns
        cwz(1,k) = -0.5*hh + hh*(k-0.5)/ns;
        cwz(2,k) = hh/ns;
        if cwz(1,k) < 0
            el(1,k) = eps1;
        else
            el(1,k) = eps2;
        end
    end
    
    MG = zeros(3,ns,N);
    for k=1:ns        
		tvar1 = ab*sign(cwz(1,k));
        tvar = 0.5*ab*hh*kg*(1 - 2*abs(cwz(1,k))/hh);
		MG(1,k,1) = 1.0;
        MG(1,k,2) = -0.5*tvar1;
        MG(1,k,N) = MG(1,k,2);
		tvar1 = tvar1/tvar; tvar2 = sqrt(1 + tvar*tvar);
        tvar3 = (tvar2 - 1.)/tvar; tvar2 = 1/tvar2;
        tvar4 = tvar3*tvar3;
		MG(2,k,1) = tvar2;
        
        tvar = tvar4;
        m = 1;
        while (2*m < no)
			MG(2,k,2*m+1) = tvar2*tvar; % 2*m
            MG(2,k,N-2*m+1) = MG(2,k,2*m+1);
            MG(3,k,2*m+1) = 1i*tvar1*tvar;
            MG(3,k,N-2*m+1) = -MG(3,k,2*m+1);
			tvar = tvar*tvar4;
            m = m + 1;
        end
        tvar = tvar3;
        % {
        m = 0;
		while ((2*m+1) < no)
			MG(2,k,2*m+2) = -tvar1*tvar;
            MG(2,k,N-2*m) = MG(2,k,2*m+2);
			MG(3,k,2*m+2) = -1i*tvar2*tvar;
            MG(3,k,N-2*m) = -MG(3,k,2*m+2);
			tvar = tvar*tvar4;
            m = m + 1;
        end
        
        MG(1,k,:) = fft(MG(1,k,:),N);
        MG(2,k,:) = fft(MG(2,k,:),N);
        MG(3,k,:) = fft(MG(3,k,:),N);
    end
end


