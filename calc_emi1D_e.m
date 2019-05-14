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
function [W] = calc_emi1D_e(V)
    global no;
    global ns;
    global eb;
    global ky;
    global kz;
    global el;
    global MG;
    global cwz;
    
    W = zeros(2*no*ns,1);
    N = 2^nextpow2(2*no-1);
    
    for k=1:ns
        kk = (k-1)*2*no;
        ff = zeros(3,no);
        ta = zeros(3,N);
        for i=1:no
            ta(1,i) = V((k-1)*2*no+i,1) + V((k-1)*2*no+i+no,1);
            ff(1,i) = -eb*ta(1,i);
            ta(3,i) = -ky(1,i)*ta(1,i);
            ff(3,i) = ta(3,i);
            ta(2,i) = kz(1,i)*( V((k-1)*2*no+i,1) - V((k-1)*2*no+i+no,1) );
            ff(2,i) = -ta(2,i);
        end
        ta = fft(ta,N,2);
        for i=1:N
            ta(1,i) = ta(1,i)*MG(1,k,i);
            tvar = ta(2,i);
            ta(2,i) = tvar*MG(2,k,i) + ta(3,i)*MG(3,k,i);
            ta(3,i) = tvar*MG(3,k,i) - ta(3,i)*MG(2,k,i);
        end
        ta = ifft(ta,N,2);
        tvar = -0.5*1i*cwz(2,k);
        for i=1:no
            ff(1,i) = ff(1,i) + el(1,k)*ta(1,i);
            ff(2,i) = ff(2,i) + ta(2,i);
            ff(3,i) = ff(3,i) + ta(3,i);
            W((k-1)*2*no+i,1) = ( ky(1,i)*ff(3,i) - ff(1,i) )/kz(1,i);
            W((k-1)*2*no+i+no,1) = tvar*( W((k-1)*2*no+i,1) + ff(2,i) );
            W((k-1)*2*no+i,1) = tvar*( W((k-1)*2*no+i,1) - ff(2,i) );
        end
    end
end

