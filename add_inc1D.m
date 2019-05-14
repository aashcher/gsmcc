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
function [SW] = add_inc1D(VI, SV)
    global no;
    global kz;
    global hh;
    global SI1;
    global SI2;
    
    SD = zeros(2,2);
    SW = zeros(2,no);
    
    for i=1:no
        texp = exp(1i*kz(1,i)*hh);
        tvar = 1/(1 - texp*texp*SI1(2,2,i)*SI2(1,1,i));
        SD(1,1) = SI1(1,1,i) + tvar*SI1(1,2,i)*SI2(1,1,i)*SI1(2,1,i)*texp*texp;
        SD(1,2) = tvar*SI1(1,2,i)*SI2(1,2,i)*texp;
        SD(2,1) = tvar*SI1(2,1,i)*SI2(2,1,i)*texp;
        SD(2,2) = SI2(2,2,i) + tvar*SI2(2,1,i)*SI1(2,2,i)*SI2(1,2,i)*texp*texp;
        
		SW(1,i) = SV(1,i) + VI(1,i)*SD(1,1) + VI(2,i)*SD(2,1);
		SW(2,i) = SV(2,i) + VI(1,i)*SD(1,2) + VI(2,i)*SD(2,2);
    end
end
