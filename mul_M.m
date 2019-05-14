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
function [W] = mul_M(V)
    global no;
    global ns;
    global MG;
    
    W = zeros(3*no*ns,1);
    N = 2^nextpow2(2*no-1);
    for k = 1:ns
        kk = (k-1)*3*no;
        ta = zeros(2,N);
        
        W((kk+1):(kk+no),1) = V((kk+1):(kk+no),1);
        ta(1,1:no) = V((kk+no+1):(kk+2*no),1);
        ta(2,1:no) = V((kk+2*no+1):(kk+3*no),1);
        
        ta = fft(ta,N,2);
        for i = 1:N
            ta(1,i) = ta(1,i)*MG(2,k,i);
            ta(2,i) = ta(2,i)*MG(2,k,i);
        end
        ta = ifft(ta,N,2);
        ta(:,no+1:N) = 0;
        ta = fft(ta,N,2);
        for i = 1:N
            ta(1,i) = ta(1,i)*MG(2,k,i);
            ta(2,i) = ta(2,i)*MG(2,k,i);
        end
        ta = ifft(ta,N,2);

        for i = 1:no
            W(kk+no+i,1) = V(kk+no+i,1) + ta(1,i);
            W(kk+2*no+i,1) = V(kk+2*no+i,1) + ta(2,i);
        end
    end
end

