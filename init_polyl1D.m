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
function init_polyl1D(kg, YZ, ab, eps1, eps2)
    global no;
    global ns;
    global hh;
    global cwz;
    global el;
    global MG;
    
    np = length(YZ(1,:));
    N = 2^nextpow2(2*no-1);
    
    z1 = max(YZ(2,:));
    z2 = min(YZ(2,:));
    kh = z1 - z2;
    hh = kh/ab;
    YZ(2,:) = (YZ(2,:) - 0.5*(z1 + z2))/kh;
    
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
    
    g = kh*kg*0.5/pi;
    MG = zeros(2,ns,N);
	for i = 1:(np-1)
		dz = YZ(2,i+1) - YZ(2,i); zz = YZ(2,i+1) + YZ(2,i);
        dy = YZ(1,i+1) - YZ(1,i); yy = YZ(1,i+1) + YZ(1,i);
		MG(1,1,1) = MG(1,1,1) + dy*zz;
		for m = 2:no
        	te = exp(-1i*pi*(m-1)*yy);
            tx = pi*(m-1)*dy;
    		tc = te*dy*(zz*sin(tx) - 1i*dz*(sin(tx)/tx - cos(tx)))/tx;
            MG(1,1,m) = MG(1,1,m) + tc;
            MG(1,1,N-m+2) = MG(1,1,N-m+2) + conj(tc);
			tc = te*dz*sin(tx)/tx;
            MG(2,1,m) = MG(2,1,m) + tc;
            MG(2,1,N-m+2) = MG(2,1,N-m+2) + conj(tc);
        end
    end
        
    for k = 2:ns
        MG(1,k,:) = MG(1,1,:);
        MG(2,k,:) = MG(2,1,:);
    end
    
    for k = 1:ns
		ts = ab*sign(cwz(1,k));
        MG(1,k,:) = -MG(1,k,:)*ts;
        ts = g*(1 - 2*abs(cwz(1,k))/hh);
        MG(2,k,:) = MG(2,k,:)*ts;
        MG(1,k,1) = MG(1,k,1) + 1;
        
        MG(1,k,:) = fft(MG(1,k,:),N);
        MG(2,k,:) = fft(MG(2,k,:),N);
    end
end
