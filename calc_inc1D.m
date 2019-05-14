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
function [V] = calc_inc1D(VI)
    global no;
    global ns;
    global hh;
    global cwz;
    global kz;
    global SI1;
    global SI2;
    
    V = zeros(2*no*ns,1);
    
    for i = 1:no
		tvar2 = 1/(1 - SI1(2,2,i)*SI2(1,1,i)*exp(2*1i*kz(1,i)*hh));
        texp = exp(1i*kz(1,i)*hh);
		tvar1 = tvar2*( VI(2,i)*SI1(2,2,i)*texp*SI2(2,1,i) + VI(1,i)*SI1(1,2,i) );
		tvar2 = tvar2*( VI(2,i)*SI2(2,1,i) + VI(1,i)*SI2(1,1,i)*texp*SI1(1,2,i) );
		for k=1:ns
			V((k-1)*2*no+i,1) = tvar1*exp(1i*kz(1,i)*(cwz(1,k)+0.5*hh));
			V((k-1)*2*no+i+no,1) = tvar2*exp(1i*kz(1,i)*(0.5*hh-cwz(1,k)));
        end
    end
end
