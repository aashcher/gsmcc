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
function [W] = calc_rad1D(V)
    global no;
    global ns;
    global hh;
    global cwz;
    global kz;
    global SI1;
    global SI2;
    
    W = zeros(2*no*ns,1);
    
    for i=1:no
        b0 = 0; b1 = 0;
        for k=1:ns
            b0 = b0 + V((k-1)*2*no+i,1) * exp(1i*kz(1,i)*(0.5*hh - cwz(1,k)));
            b1 = b1 + V((k-1)*2*no+i+no,1) * exp(1i*kz(1,i)*(0.5*hh + cwz(1,k)));
        end
		tvar1 = exp(1i*kz(1,i)*hh);
		tvar2 = b0;
        b0 = b1 + b0*SI2(1,1,i)*tvar1;
        b1 = tvar2 + b1*SI1(2,2,i)*tvar1;
		tvar2 = 1/(1 - tvar1*tvar1*SI1(2,2,i)*SI2(1,1,i));
        b0 = b0*SI1(2,2,i)*tvar2*exp(1i*kz(1,i)*0.5*cwz(2,1));
        b1 = b1*SI2(1,1,i)*tvar2*exp(1i*kz(1,i)*0.5*cwz(2,ns));
		for k=1:(ns-1)
			tvar1 = V((k-1)*2*no+i,1);
            W((k-1)*2*no+i,1) = 0.5*tvar1 + b0;
			b0 = (b0 + tvar1)*exp(1i*kz(1,i)*cwz(2,k));
			tvar1 = V((ns-k)*2*no+i+no,1);
            W((ns-k)*2*no+i+no,1) = 0.5*tvar1 + b1;
			b1 = (b1 + tvar1)*exp(1i*kz(1,i)*cwz(2,ns-k+1));
        end
		W((ns-1)*2*no+i,1) = 0.5*V((ns-1)*2*no+i,1) + b0;
		W(i+no,1) = 0.5*V(i+no,1) + b1;
    end
end

